"""
Decode CSSHoneycomb memory experiments with PyMatching and estimate a logical error rate.

Experiment
----------
    prepare -> bZ -> rX -> gZ -> bX -> (perfect final readout of the whole ISG)

The final readout (`CSSHoneycomb.final_readout`) is noiseless: it measures every
plaquette stabilizer that is still open with MPP, closing the open detectors as if
the circuit were going to carry on into further rounds.

Typical use
-----------
    from decoding import run_experiment
    result = run_experiment(config_filepath="config/5x3.json", max_shots=200_000)
    print(result)

or, piece by piece:

    circuit  = build_circuit(config_filepath="config/5x3.json")
    check_circuit(circuit)
    matching = build_matching(circuit)
    result   = estimate_logical_error_rate(circuit, matching, max_shots=200_000)
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field

import numpy as np
import pymatching
import stim

COLOUR_MAP = {"r": "red", "g": "green", "b": "blue"}

DEFAULT_CYCLE = ("bZ", "rX", "gZ", "bX")


# --------------------------------------------------------------------------- #
# circuit construction
# --------------------------------------------------------------------------- #

def measurement_rounds(honey, cycle) -> str:
    """Measure each (colour, flavour) round in `cycle`, in order. Mirrors the notebook helper."""
    string = ""
    for token in cycle:
        colour = COLOUR_MAP[token[0]]
        flavour = token[-1]

        string += f"#{colour}{flavour}\n"
        string += honey.get_measurement(colour, flavour) + "\n"
        string += honey.tick()
    return string


def build_circuit(
    cycle=DEFAULT_CYCLE,
    *,
    ncols: int = 5,
    nrows: int = 3,
    config_filepath: str = "config/5x3.json",
    error_model: str | None = None,
    bases: tuple[str, ...] = ("Z", "X"),
    final_readout: bool = True,
) -> stim.Circuit:
    """
    Build the noisy memory circuit.

    Parameters
    ----------
    cycle : sequence of str
        Rounds to measure, e.g. ("bZ", "rX", "gZ", "bX"). The logical observables are
        evaluated after the LAST round in this list.
    ncols, nrows, config_filepath, error_model
        Passed straight to `CSSHoneycomb`.
    bases : tuple[str, ...]
        Logicals to track. Observable i corresponds to bases[i]; the default gives
        OBSERVABLE 0 = Z logical, OBSERVABLE 1 = X logical.
    final_readout : bool
        Append the noiseless final readout of all remaining ISG stabilizers.

    Notes
    -----
    A fresh `CSSHoneycomb` is created on every call because it is stateful
    (clock, measurement list, detector bookkeeping), so it cannot be reused.
    """
    from honeycomb import CSSHoneycomb  # lazy: keeps this module importable without the lattice code

    cycle = list(cycle)

    s = CSSHoneycomb(ncols=ncols, nrows=nrows, error_model=error_model,
                     config_filepath=config_filepath)
    s.set_measurement_cycle(cycle)

    string = s.initialise_circuit()
    string += s.prepare_qubits()
    string += measurement_rounds(s, cycle)

    # The observable's flow is solved on a noiseless copy: the measurement records it
    # finds are the same once noise is added, and this avoids stim's flow solver ever
    # seeing noise channels.
    last_colour = COLOUR_MAP[cycle[-1][0]]
    last_flavour = cycle[-1][-1]
    noiseless = str(stim.Circuit(string).without_noise())
    string += s.write_logical_observable(noiseless, last_colour, last_flavour, bases=bases)

    if final_readout:
        string += s.final_readout()

    return stim.Circuit(string)


def check_circuit(circuit: stim.Circuit, shots: int = 64) -> None:
    """
    Sanity-check the circuit before decoding.

    Samples the noiseless circuit and verifies that every detector and observable is
    deterministically 0. If this fails, no decoder will give a meaningful answer:
    either a detector is wrong, or an observable has a constant offset.
    """
    dets, obs = circuit.without_noise().compile_detector_sampler().sample(
        shots, separate_observables=True
    )
    bad_dets = np.flatnonzero(dets.any(axis=0))
    bad_obs = np.flatnonzero(obs.any(axis=0))
    if bad_dets.size or bad_obs.size:
        raise ValueError(
            "Noiseless circuit is not deterministic.\n"
            f"  detectors that fire: {bad_dets.tolist()}\n"
            f"  observables that flip: {bad_obs.tolist()}"
        )


# --------------------------------------------------------------------------- #
# decoder construction
# --------------------------------------------------------------------------- #

def _split_components(targets):
    """Split DEM error targets on `^` into the graphlike components."""
    components = [[]]
    for t in targets:
        if t.is_separator():
            components.append([])
        else:
            components[-1].append(t)
    return components


def _drop_hyperedges(dem: stim.DetectorErrorModel):
    """Remove error mechanisms that still touch more than two detectors."""
    out = stim.DetectorErrorModel()
    dropped = 0
    dropped_p = 0.0
    for inst in dem.flattened():
        if inst.type == "error":
            components = _split_components(inst.targets_copy())
            if any(sum(t.is_relative_detector_id() for t in c) > 2 for c in components):
                dropped += 1
                dropped_p += inst.args_copy()[0]
                continue
        out.append(inst)
    return out, dropped, dropped_p


def build_dem(circuit: stim.Circuit, *, drop_hyperedges: bool = True) -> stim.DetectorErrorModel:
    """
    Detector error model with graphlike (<= 2 detector) error mechanisms.

    First tries stim's own decomposition. If some errors cannot be decomposed, and
    `drop_hyperedges` is True, those mechanisms are removed from the DEM the decoder
    sees (a warning reports how many). The circuit that is *sampled* is unchanged, so
    the logical error rate you measure is still honest for this decoder; it is just a
    weaker decoder than one that could use the hyperedges.
    """
    kwargs = dict(decompose_errors=True, approximate_disjoint_errors=True)
    try:
        return circuit.detector_error_model(**kwargs)
    except ValueError as exc:
        if "decompose" not in str(exc).lower() or not drop_hyperedges:
            raise
    dem = circuit.detector_error_model(ignore_decomposition_failures=True, **kwargs)
    dem, dropped, dropped_p = _drop_hyperedges(dem)
    warnings.warn(
        f"{dropped} error mechanisms could not be decomposed into graphlike pieces and were "
        f"removed from the matching graph (total probability {dropped_p:.3g}).",
        stacklevel=2,
    )
    return dem


def _merge_undetectable_components(dem: stim.DetectorErrorModel) -> stim.DetectorErrorModel:
    """
    Correlated matching rejects decomposed errors with a component that has no detectors
    (e.g. `D0 D1 ^ L0`). Fold the observables of such components into a detected component
    of the same error, so the error keeps its total effect. Errors that are entirely
    undetectable are collapsed to a single component.

    This is a mild approximation: the folded observable is attributed to one matching edge.
    It only affects the few mechanisms that flip a logical alongside a detected fault.
    """
    def has_det(c):
        return any(t.is_relative_detector_id() for t in c)

    def obs_ids(c):
        return {t.val for t in c if t.is_logical_observable_id()}

    out = stim.DetectorErrorModel()
    for inst in dem.flattened():
        if inst.type != "error":
            out.append(inst)
            continue
        comps = _split_components(inst.targets_copy())
        if len(comps) == 1 or all(has_det(c) for c in comps):
            out.append(inst)
            continue

        detected = [c for c in comps if has_det(c)]
        extra_obs = set()
        for c in comps:
            if not has_det(c):
                extra_obs ^= obs_ids(c)

        targets = []
        if not detected:                      # fully undetectable: one component
            for i in sorted(extra_obs):
                targets.append(stim.target_logical_observable_id(i))
        else:
            for k, c in enumerate(detected):
                if k:
                    targets.append(stim.target_separator())
                dets = [t for t in c if t.is_relative_detector_id()]
                obs = obs_ids(c) ^ (extra_obs if k == 0 else set())
                targets += [stim.target_relative_detector_id(t.val) for t in dets]
                targets += [stim.target_logical_observable_id(i) for i in sorted(obs)]
        out.append("error", inst.args_copy(), targets)
    return out


def _canonicalise_edges(dem: stim.DetectorErrorModel) -> stim.DetectorErrorModel:
    """
    Merge parallel matching edges ourselves, keeping the most probable observable class.

    PyMatching merges parallel edges but keeps the observables of whichever edge it saw FIRST,
    not the most likely one, so the decoder's behaviour depends on the order of the DEM. Here every
    graphlike component becomes one edge per detector set: its probability is the independent
    combination of all classes on that edge, and its observables are those of the class with the
    largest probability. Components with no detectors cannot be matched and are skipped.

    This discards the X/Z correlation information, so it is for plain (uncorrelated) matching.
    """
    from collections import defaultdict

    classes = defaultdict(lambda: defaultdict(float))     # detector set -> observable set -> p
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        p = inst.args_copy()[0]
        for comp in _split_components(inst.targets_copy()):
            dets = frozenset(t.val for t in comp if t.is_relative_detector_id())
            if not dets:
                continue
            obs = frozenset(t.val for t in comp if t.is_logical_observable_id())
            q = classes[dets][obs]
            classes[dets][obs] = q * (1 - p) + p * (1 - q)

    out = stim.DetectorErrorModel()
    for dets, by_obs in classes.items():
        total = 0.0
        for q in by_obs.values():
            total = total * (1 - q) + q * (1 - total)
        best = max(by_obs, key=by_obs.get)
        targets = [stim.target_relative_detector_id(d) for d in sorted(dets)]
        targets += [stim.target_logical_observable_id(o) for o in sorted(best)]
        out.append("error", total, targets)
    for i in range(dem.num_detectors):
        out.append("detector", [], [stim.target_relative_detector_id(i)])
    for i in range(dem.num_observables):
        out.append("logical_observable", [], [stim.target_logical_observable_id(i)])
    return out


def build_matching(
    circuit: stim.Circuit, *, drop_hyperedges: bool = True, correlated: bool = False
) -> pymatching.Matching:
    """
    PyMatching decoder for the circuit.

    correlated=True builds the matcher for correlated matching (needs pymatching >= 2.3),
    which uses the X/Z decomposition of Y-type errors instead of decoding the two halves
    independently. Pass the same flag to `estimate_logical_error_rate`.
    """
    dem = build_dem(circuit, drop_hyperedges=drop_hyperedges)
    if correlated:
        dem = _merge_undetectable_components(dem)
        return pymatching.Matching.from_detector_error_model(dem, enable_correlations=True)
    return pymatching.Matching.from_detector_error_model(_canonicalise_edges(dem))


def first_order_failure_rate(circuit: stim.Circuit, matching: pymatching.Matching,
                             correlated: bool = False, top: int = 0):
    """
    Exact first-order logical failure rate of THIS decoder: decode every single-fault
    syndrome in the (undecomposed) DEM and sum the probabilities of the ones it gets wrong.
    Compare with the ideal-decoder estimate to see how much the decoder itself costs.

    Returns (rate, worst) where `worst` is the `top` failing mechanisms as (p, DEM text).
    """
    dem = circuit.detector_error_model(approximate_disjoint_errors=True)
    probs, rows = [], []
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        ts = inst.targets_copy()
        rows.append((
            [t.val for t in ts if t.is_relative_detector_id()],
            [t.val for t in ts if t.is_logical_observable_id()],
        ))
        probs.append(inst.args_copy()[0])

    dets = np.zeros((len(rows), circuit.num_detectors), dtype=bool)
    obs = np.zeros((len(rows), circuit.num_observables), dtype=bool)
    for i, (d, o) in enumerate(rows):
        dets[i, d] = True
        obs[i, o] = True

    kw = {"enable_correlations": True} if correlated else {}
    pred = np.asarray(matching.decode_batch(dets, **kw), dtype=bool)
    if pred.shape[1] < obs.shape[1]:
        pred = np.pad(pred, ((0, 0), (0, obs.shape[1] - pred.shape[1])))
    fail = (pred[:, :obs.shape[1]] != obs).any(axis=1)

    probs = np.asarray(probs)
    worst = []
    if top:
        for i in np.argsort(-(probs * fail))[:top]:
            if fail[i]:
                worst.append((float(probs[i]), f"D{rows[i][0]} L{rows[i][1]}"))
    return float(probs[fail].sum()), worst


# --------------------------------------------------------------------------- #
# sampling + decoding
# --------------------------------------------------------------------------- #

def wilson_interval(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion (default 95%)."""
    if n == 0:
        return 0.0, 1.0
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


@dataclass
class DecodingResult:
    shots: int
    failures: int                                # shots where ANY tracked logical was wrong
    failures_per_observable: list[int]
    observable_labels: tuple[str, ...] = field(default=())

    @property
    def logical_error_rate(self) -> float:
        return self.failures / self.shots if self.shots else float("nan")

    @property
    def confidence_interval(self) -> tuple[float, float]:
        return wilson_interval(self.failures, self.shots)

    def rate_for(self, label: str) -> float:
        i = self.observable_labels.index(label)
        return self.failures_per_observable[i] / self.shots

    def __str__(self) -> str:
        lo, hi = self.confidence_interval
        lines = [
            f"shots: {self.shots}",
            f"logical error rate (any logical): {self.logical_error_rate:.4g}"
            f"  [95% CI {lo:.4g}, {hi:.4g}]  ({self.failures} failures)",
        ]
        for label, k in zip(self.observable_labels, self.failures_per_observable):
            lines.append(f"  {label} logical: {k / self.shots:.4g}  ({k} failures)")
        return "\n".join(lines)


def estimate_logical_error_rate(
    circuit: stim.Circuit,
    matching: pymatching.Matching | None = None,
    *,
    max_shots: int = 100_000,
    max_errors: int | None = None,
    batch_size: int = 10_000,
    seed: int | None = None,
    observable_labels: tuple[str, ...] = ("Z", "X"),
    correlated: bool = False,
) -> DecodingResult:
    """
    Sample the circuit, decode with PyMatching, and count logical failures.

    A shot fails if the decoder's prediction for ANY tracked observable differs from
    the sampled value. Per-observable counts are reported as well.

    Parameters
    ----------
    max_shots : int
        Upper limit on the number of shots.
    max_errors : int | None
        Stop early once this many failing shots have been seen (useful at high error
        rates, where you don't need all of `max_shots`).
    batch_size : int
        Shots sampled and decoded per batch.
    observable_labels : tuple[str, ...]
        Names for the observables, in order; should match the `bases` used to build the circuit.
    """
    if matching is None:
        matching = build_matching(circuit, correlated=correlated)
    decode_kw = {"enable_correlations": True} if correlated else {}

    num_obs = circuit.num_observables
    if len(observable_labels) != num_obs:
        observable_labels = tuple(f"obs{i}" for i in range(num_obs))

    sampler = circuit.compile_detector_sampler(seed=seed)

    shots = 0
    failures = 0
    per_obs = np.zeros(num_obs, dtype=np.int64)

    while shots < max_shots and (max_errors is None or failures < max_errors):
        n = min(batch_size, max_shots - shots)
        dets, obs = sampler.sample(n, separate_observables=True)

        pred = np.asarray(matching.decode_batch(dets, **decode_kw), dtype=bool)
        if pred.ndim == 1:
            pred = pred[:, None]
        if pred.shape[1] < num_obs:      # observables that never appear in the DEM
            pred = np.pad(pred, ((0, 0), (0, num_obs - pred.shape[1])))

        wrong = pred[:, :num_obs] != obs
        failures += int(wrong.any(axis=1).sum())
        per_obs += wrong.sum(axis=0)
        shots += n

    return DecodingResult(
        shots=shots,
        failures=failures,
        failures_per_observable=per_obs.tolist(),
        observable_labels=tuple(observable_labels),
    )


# --------------------------------------------------------------------------- #
# one-call convenience
# --------------------------------------------------------------------------- #

def run_experiment(
    cycle=DEFAULT_CYCLE,
    *,
    ncols: int = 5,
    nrows: int = 3,
    config_filepath: str = "config/5x3.json",
    error_model: str | None = None,
    bases: tuple[str, ...] = ("Z", "X"),
    max_shots: int = 100_000,
    max_errors: int | None = None,
    batch_size: int = 10_000,
    seed: int | None = None,
    check: bool = True,
    correlated: bool = False,
) -> DecodingResult:
    """Build the circuit, (optionally) sanity-check it, decode, and return the result."""
    circuit = build_circuit(
        cycle, ncols=ncols, nrows=nrows, config_filepath=config_filepath,
        error_model=error_model, bases=bases,
    )
    if check:
        check_circuit(circuit)
    return estimate_logical_error_rate(
        circuit, build_matching(circuit, correlated=correlated),
        max_shots=max_shots, max_errors=max_errors, batch_size=batch_size,
        seed=seed, observable_labels=tuple(bases), correlated=correlated,
    )


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config/5x3.json")
    ap.add_argument("--error-model", default=None)
    ap.add_argument("--ncols", type=int, default=5)
    ap.add_argument("--nrows", type=int, default=3)
    ap.add_argument("--cycle", nargs="+", default=list(DEFAULT_CYCLE))
    ap.add_argument("--shots", type=int, default=100_000)
    ap.add_argument("--max-errors", type=int, default=None)
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    print(run_experiment(
        args.cycle, ncols=args.ncols, nrows=args.nrows, config_filepath=args.config,
        error_model=args.error_model, max_shots=args.shots, max_errors=args.max_errors,
        seed=args.seed,
    ))