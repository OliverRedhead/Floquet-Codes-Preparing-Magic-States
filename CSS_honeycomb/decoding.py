"""
Decode CSSHoneycomb memory experiments and estimate logical error rates, using sinter.

Experiment
----------
    prepare -> bZ -> rX -> gZ -> bX -> (perfect final readout of the whole ISG)

The final readout (`CSSHoneycomb.final_readout`) is noiseless: it measures every
plaquette stabilizer that is still open with MPP, closing the open detectors as if
the circuit were going to carry on into further rounds.

Decoders
--------
    "pymatching"             plain matching on a graphlike DEM (hyperedges decomposed or dropped,
                             parallel edges merged deterministically)
    "pymatching-correlated"  correlated matching (pymatching >= 2.3)
    "bposd"                  BP+OSD on the undecomposed DEM (needs `pip install stimbposd`)

Each decoder is given the DEM it needs: the matching decoders get a graphlike DEM, BP+OSD gets
the full hyperedge DEM. Comparing decoders on the same circuit is the cleanest way to tell a
physical effect from a decomposition artifact.

Typical use
-----------
    from decoding import run_experiment, format_results

    results = run_experiment(
        config_filepath="config/5x3.json",
        decoders=("pymatching", "pymatching-correlated", "bposd"),
        max_shots=1_000_000, max_errors=2000, num_workers=8,
    )
    print(format_results(results))

or, to sweep several circuits in one sinter run (with resume support):

    tasks  = make_tasks(circuit_a, metadata={"rep": "A"}) + make_tasks(circuit_b, metadata={"rep": "B"})
    stats  = collect(tasks, max_shots=1_000_000, max_errors=2000, save_resume_filepath="out.csv")
    print(format_results(to_results(stats)))

On macOS/Windows, call sinter from under `if __name__ == "__main__":` (multiprocessing uses spawn).
"""

from __future__ import annotations

import math
import os
import warnings
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
import pymatching
import sinter
import stim

COLOUR_MAP = {"r": "red", "g": "green", "b": "blue"}

DEFAULT_CYCLE = ("bZ", "rX", "gZ", "bX")

MATCHING = "pymatching"
CORRELATED = "pymatching-correlated"
BPOSD = "bposd"


# --------------------------------------------------------------------------- #
# circuit construction
# --------------------------------------------------------------------------- #

def measurement_rounds(honey, cycle) -> str:
    """Measure each (colour, flavour) round in `cycle`, in order."""
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

    Observable i corresponds to bases[i]; the default gives OBSERVABLE 0 = Z logical,
    OBSERVABLE 1 = X logical. The observables are evaluated after the LAST round in `cycle`.

    A fresh `CSSHoneycomb` is created on every call because it is stateful (clock,
    measurement list, detector bookkeeping), so it cannot be reused.
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
    # finds are the same once noise is added, and this keeps stim's flow solver away
    # from noise channels.
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

    1. Asks stim to build a DEM of the noiseless circuit, which raises if any detector or
       observable is non-deterministic (anti-commuting with the stabilizer group).
    2. Samples the noiseless circuit and verifies every detector and observable is 0, which
       catches constant offsets that the DEM check cannot see.
    """
    if circuit.num_observables == 0:
        raise ValueError("Circuit has no OBSERVABLE_INCLUDE instructions.")

    noiseless = circuit.without_noise()
    try:
        noiseless.detector_error_model()
    except ValueError as exc:
        raise ValueError(f"A detector or observable is non-deterministic: {exc}") from exc

    dets, obs = noiseless.compile_detector_sampler().sample(shots, separate_observables=True)
    bad_dets = np.flatnonzero(dets.any(axis=0))
    bad_obs = np.flatnonzero(obs.any(axis=0))
    if bad_dets.size or bad_obs.size:
        raise ValueError(
            "Noiseless circuit has deterministic non-zero outcomes.\n"
            f"  detectors that fire: {bad_dets.tolist()}\n"
            f"  observables that flip: {bad_obs.tolist()}"
        )


def single_observable_circuit(circuit: stim.Circuit, index: int) -> stim.Circuit:
    """Copy of `circuit` that keeps only OBSERVABLE_INCLUDE(index), renumbered to 0."""
    out = stim.Circuit()
    for inst in circuit.flattened():
        if inst.name == "OBSERVABLE_INCLUDE":
            if int(inst.gate_args_copy()[0]) != index:
                continue
            out.append("OBSERVABLE_INCLUDE", inst.targets_copy(), 0)
        else:
            out.append(inst)
    return out


# --------------------------------------------------------------------------- #
# detector error models
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

    First tries stim's own decomposition. If some errors cannot be decomposed and
    `drop_hyperedges` is True, those mechanisms are removed from the DEM the decoder sees
    (a warning reports how many). The circuit that is *sampled* is unchanged, so the logical
    error rate is still honest for this decoder; it is just weaker than one using hyperedges.
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
    of the same error so the error keeps its total effect. Fully undetectable errors are
    collapsed to a single component.

    This is a mild approximation: the folded observable is attributed to one matching edge.
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
    so its behaviour depends on the order of the DEM (and therefore on which logical
    representative you chose). Here every graphlike component becomes one edge per detector
    set: its probability is the independent combination of all classes on that edge, and its
    observables are those of the class with the largest probability. Components with no
    detectors cannot be matched and are skipped.

    This discards X/Z correlation information, so it is for plain (uncorrelated) matching.
    """
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


def decoder_dem(
    circuit: stim.Circuit,
    decoder: str,
    *,
    drop_hyperedges: bool = True,
    _graphlike: stim.DetectorErrorModel | None = None,
) -> stim.DetectorErrorModel:
    """The DEM that `decoder` should be compiled from."""
    if decoder in (MATCHING, CORRELATED):
        graphlike = _graphlike if _graphlike is not None else build_dem(circuit, drop_hyperedges=drop_hyperedges)
        return _canonicalise_edges(graphlike) if decoder == MATCHING else _merge_undetectable_components(graphlike)
    # everything else (e.g. BP+OSD) works on the full, undecomposed DEM
    return circuit.detector_error_model(approximate_disjoint_errors=True)


def build_matching(
    circuit: stim.Circuit, *, drop_hyperedges: bool = True, correlated: bool = False
) -> pymatching.Matching:
    """
    In-process PyMatching decoder for the circuit (used by `estimate_logical_error_rate` and
    `first_order_failure_rate`). For large runs use `make_tasks` + `collect` instead.
    """
    dem = decoder_dem(circuit, CORRELATED if correlated else MATCHING, drop_hyperedges=drop_hyperedges)
    if correlated:
        return pymatching.Matching.from_detector_error_model(dem, enable_correlations=True)
    return pymatching.Matching.from_detector_error_model(dem)


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
    if pred.ndim == 1:
        pred = pred[:, None]
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
# sinter decoders
# --------------------------------------------------------------------------- #

class _CompiledCorrelatedMatching(sinter.CompiledDecoder):
    def __init__(self, matching: pymatching.Matching, num_detectors: int, num_observables: int):
        self.matching = matching
        self.num_detectors = num_detectors
        self.num_observables = num_observables

    def decode_shots_bit_packed(self, *, bit_packed_detection_event_data: np.ndarray) -> np.ndarray:
        dets = np.unpackbits(
            bit_packed_detection_event_data, axis=1, bitorder="little", count=self.num_detectors
        ).astype(bool)
        pred = np.asarray(self.matching.decode_batch(dets, enable_correlations=True), dtype=np.uint8)
        if pred.ndim == 1:
            pred = pred[:, None]
        if pred.shape[1] < self.num_observables:      # observables that never appear in the DEM
            pred = np.pad(pred, ((0, 0), (0, self.num_observables - pred.shape[1])))
        return np.packbits(pred[:, :self.num_observables], axis=1, bitorder="little")


class CorrelatedMatchingDecoder(sinter.Decoder):
    """Sinter wrapper for PyMatching's correlated matching (needs pymatching >= 2.3)."""

    def compile_decoder_for_dem(self, *, dem: stim.DetectorErrorModel) -> sinter.CompiledDecoder:
        matching = pymatching.Matching.from_detector_error_model(dem, enable_correlations=True)
        return _CompiledCorrelatedMatching(matching, dem.num_detectors, dem.num_observables)


def custom_decoders() -> dict[str, sinter.Decoder]:
    """Decoders sinter doesn't ship with. BP+OSD is included only if `stimbposd` is installed."""
    decoders: dict[str, sinter.Decoder] = {CORRELATED: CorrelatedMatchingDecoder()}
    try:
        from stimbposd import sinter_decoders
        decoders.update(sinter_decoders())
    except ImportError:
        pass
    return decoders


# --------------------------------------------------------------------------- #
# sinter tasks, collection, results
# --------------------------------------------------------------------------- #

def make_tasks(
    circuit: stim.Circuit,
    decoders: tuple[str, ...] = (MATCHING,),
    *,
    observable_labels: tuple[str, ...] = ("Z", "X"),
    per_observable: bool = False,
    drop_hyperedges: bool = True,
    metadata: dict | None = None,
) -> list[sinter.Task]:
    """
    One sinter Task per decoder, each carrying the DEM that decoder should use.

    A shot counts as a failure if ANY observable is mispredicted (metadata "observable": "any").

    per_observable=True adds tasks that track a single logical at a time (metadata "observable":
    "Z"/"X"). Those decode a DEM containing only that observable, so they are a slightly
    different decoder from the joint one (matching picks edge observables per observable
    rather than jointly) and cost extra shots. For a breakdown of the joint decoder's
    failures use `estimate_logical_error_rate` instead.

    `metadata` is merged into each task's json_metadata, so you can label circuits
    (e.g. {"rep": "A", "p": 1e-3}) and collect several circuits in one `collect` call.
    """
    metadata = dict(metadata or {})
    labels = tuple(observable_labels)
    if len(labels) != circuit.num_observables:
        labels = tuple(f"obs{i}" for i in range(circuit.num_observables))

    variants = [("any", circuit)]
    if per_observable and circuit.num_observables > 1:
        variants += [(lab, single_observable_circuit(circuit, k)) for k, lab in enumerate(labels)]

    tasks = []
    for obs_label, circ in variants:
        graphlike = None
        for decoder in decoders:
            if decoder in (MATCHING, CORRELATED) and graphlike is None:
                graphlike = build_dem(circ, drop_hyperedges=drop_hyperedges)
            dem = decoder_dem(circ, decoder, drop_hyperedges=drop_hyperedges, _graphlike=graphlike)
            tasks.append(sinter.Task(
                circuit=circ,
                detector_error_model=dem,
                decoder=decoder,
                json_metadata={**metadata, "observable": obs_label},
            ))
    return tasks


def collect(
    tasks: list[sinter.Task],
    *,
    num_workers: int | None = None,
    max_shots: int | None = 100_000,
    max_errors: int | None = None,
    save_resume_filepath: str | None = None,
    print_progress: bool = False,
    extra_decoders: dict[str, sinter.Decoder] | None = None,
) -> list[sinter.TaskStats]:
    """
    Run tasks in parallel with sinter.

    Each task stops at `max_shots` or `max_errors`, whichever comes first (use max_errors to
    fix the relative statistical error at about 1/sqrt(max_errors)). With
    `save_resume_filepath` the stats are appended to a CSV and tasks already finished there are
    not re-run, so interrupted sweeps can be resumed.
    """
    decoders = custom_decoders()
    decoders.update(extra_decoders or {})
    return sinter.collect(
        num_workers=num_workers or os.cpu_count() or 1,
        tasks=tasks,
        custom_decoders=decoders,
        max_shots=max_shots,
        max_errors=max_errors,
        save_resume_filepath=save_resume_filepath,
        print_progress=print_progress,
    )


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
    decoder: str
    observable: str                    # "any", or a single observable label
    shots: int
    failures: int
    discards: int = 0
    seconds: float = float("nan")
    metadata: dict = field(default_factory=dict)
    failures_per_observable: dict = field(default_factory=dict)   # only from the in-process estimator

    @property
    def logical_error_rate(self) -> float:
        return self.failures / self.shots if self.shots else float("nan")

    @property
    def confidence_interval(self) -> tuple[float, float]:
        return wilson_interval(self.failures, self.shots)

    def __str__(self) -> str:
        lo, hi = self.confidence_interval
        extra = {k: v for k, v in self.metadata.items() if k != "observable"}
        tag = f" {extra}" if extra else ""
        s = (f"{self.decoder:<22} {self.observable:<5} {self.logical_error_rate:.4g}"
             f"  [95% CI {lo:.4g}, {hi:.4g}]  ({self.failures}/{self.shots}){tag}")
        for label, k in self.failures_per_observable.items():
            s += f"\n    {label} logical: {k / self.shots:.4g}  ({k} failures)"
        return s


def to_results(stats: list[sinter.TaskStats]) -> list[DecodingResult]:
    """Convert sinter TaskStats into DecodingResults."""
    return [
        DecodingResult(
            decoder=s.decoder,
            observable=s.json_metadata.get("observable", "any"),
            shots=s.shots,
            failures=s.errors,
            discards=s.discards,
            seconds=s.seconds,
            metadata=dict(s.json_metadata),
        )
        for s in stats
    ]


def format_results(results: list[DecodingResult]) -> str:
    return "\n".join(str(r) for r in results)


# --------------------------------------------------------------------------- #
# in-process estimator (reproducible, gives per-observable counts of the joint decoder)
# --------------------------------------------------------------------------- #

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
    Single-process sample-and-decode with PyMatching. A shot fails if the prediction for ANY
    tracked observable is wrong; per-observable counts of that same decoder are also reported.
    Use `collect` for anything big: it is parallel, resumable and supports other decoders.
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
        if pred.shape[1] < num_obs:
            pred = np.pad(pred, ((0, 0), (0, num_obs - pred.shape[1])))

        wrong = pred[:, :num_obs] != obs
        failures += int(wrong.any(axis=1).sum())
        per_obs += wrong.sum(axis=0)
        shots += n

    return DecodingResult(
        decoder=CORRELATED if correlated else MATCHING,
        observable="any",
        shots=shots,
        failures=failures,
        failures_per_observable=dict(zip(observable_labels, per_obs.tolist())),
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
    decoders: tuple[str, ...] = (MATCHING,),
    per_observable: bool = False,
    max_shots: int | None = 100_000,
    max_errors: int | None = None,
    num_workers: int | None = None,
    save_resume_filepath: str | None = None,
    print_progress: bool = False,
    check: bool = True,
    metadata: dict | None = None,
) -> list[DecodingResult]:
    """Build the circuit, (optionally) sanity-check it, decode with each decoder via sinter."""
    circuit = build_circuit(
        cycle, ncols=ncols, nrows=nrows, config_filepath=config_filepath,
        error_model=error_model, bases=bases,
    )
    if check:
        check_circuit(circuit)
    tasks = make_tasks(
        circuit, decoders, observable_labels=tuple(bases),
        per_observable=per_observable, metadata=metadata,
    )
    stats = collect(
        tasks, num_workers=num_workers, max_shots=max_shots, max_errors=max_errors,
        save_resume_filepath=save_resume_filepath, print_progress=print_progress,
    )
    return to_results(stats)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config/5x3.json")
    ap.add_argument("--error-model", default=None)
    ap.add_argument("--ncols", type=int, default=5)
    ap.add_argument("--nrows", type=int, default=3)
    ap.add_argument("--cycle", nargs="+", default=list(DEFAULT_CYCLE))
    ap.add_argument("--decoders", nargs="+", default=[MATCHING], help=f"{MATCHING} {CORRELATED} {BPOSD}")
    ap.add_argument("--shots", type=int, default=100_000)
    ap.add_argument("--max-errors", type=int, default=None)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--per-observable", action="store_true")
    ap.add_argument("--save-resume", default=None, help="CSV file for sinter stats; reruns resume from it")
    args = ap.parse_args()

    print(format_results(run_experiment(
        args.cycle, ncols=args.ncols, nrows=args.nrows, config_filepath=args.config,
        error_model=args.error_model, decoders=tuple(args.decoders),
        per_observable=args.per_observable, max_shots=args.shots, max_errors=args.max_errors,
        num_workers=args.workers, save_resume_filepath=args.save_resume, print_progress=True,
    )))