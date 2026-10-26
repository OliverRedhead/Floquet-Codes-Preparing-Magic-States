"""Graphlike distance after an ideal injection into the codespace.

Four noiseless rounds (bZ, rX, gZ, bX) prepare the 5x3 lattice and inject the
T state. One six-round cycle then runs with Hadamard, CNOT, measurement, and
idle noise from error_model.json. Preparation noise stays off, including the
ancilla resets inside each measurement round, which share that channel.

The distance is the same search as testing/simulation.ipynb: drop detector-free
faults, then len(shortest_graphlike_error()).
"""

from __future__ import annotations

import sys
from pathlib import Path

import stim

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from honeycomb import CSSHoneycomb

IDEAL = ("bZ", "rX", "gZ", "bX")
NOISY = ("rZ", "gX", "bZ", "rX", "gZ", "bX")
COLOUR = {"r": "red", "g": "green", "b": "blue"}
ALL_FLAGS = (
    "initialisation_error",
    "hadamard_error",
    "cnot_error",
    "measurement_error",
    "idling_error",
)
# Ancilla resets use initialisation_error, so that flag stays off.
ROUND_FLAGS = (
    "hadamard_error",
    "cnot_error",
    "measurement_error",
    "idling_error",
)


def _set_flags(model, enabled: tuple[str, ...]) -> None:
    for key in ALL_FLAGS:
        model.config[key] = key in enabled


def _measurement_rounds(honey: CSSHoneycomb, cycle: tuple[str, ...]) -> str:
    string = ""
    for token in cycle:
        colour = COLOUR[token[0]]
        flavour = token[-1]
        string += f"#{colour}{flavour}\n"
        string += honey.get_measurement(colour, flavour) + "\n"
        string += honey.tick()
    return string


def build_circuit() -> stim.Circuit:
    """Ideal injection, then one noisy six-round cycle, then a noiseless readout."""
    honey = CSSHoneycomb(
        ncols=5,
        nrows=3,
        error_model=str(ROOT / "error_model.json"),
        config_filepath=str(ROOT / "config" / "5x3.json"),
    )
    honey.set_measurement_cycle([*IDEAL, *NOISY])

    _set_flags(honey.error_model, ())
    string = honey.initialise_circuit()
    string += honey.prepare_qubits()
    string += _measurement_rounds(honey, IDEAL)

    _set_flags(honey.error_model, ROUND_FLAGS)
    string += _measurement_rounds(honey, NOISY)

    last = NOISY[-1]
    noiseless = str(stim.Circuit(string).without_noise())
    string += honey.write_logical_observable(
        noiseless, COLOUR[last[0]], last[-1], bases=("Z", "X")
    )
    string += honey.final_readout()
    return stim.Circuit(string)


def check_noiseless(circuit: stim.Circuit, shots: int = 64) -> None:
    """Every detector and both observables are 0 when the noise is removed."""
    dets, obs = circuit.without_noise().compile_detector_sampler().sample(
        shots, separate_observables=True
    )
    bad_dets = [i for i in range(dets.shape[1]) if dets[:, i].any()]
    bad_obs = [i for i in range(obs.shape[1]) if obs[:, i].any()]
    if bad_dets or bad_obs:
        raise ValueError(
            "Noiseless circuit is not deterministic.\n"
            f"  detectors that fire: {bad_dets}\n"
            f"  observables that flip: {bad_obs}"
        )


def effective_distance(circuit: stim.Circuit) -> int:
    """Notebook search: shortest graphlike logical that still touches a detector."""
    dem = circuit.detector_error_model(
        decompose_errors=True,
        approximate_disjoint_errors=True,
        ignore_decomposition_failures=True,
    )
    clean = stim.DetectorErrorModel()
    for inst in dem.flattened():
        if inst.type == "error" and not any(
            t.is_relative_detector_id() for t in inst.targets_copy()
        ):
            continue
        clean.append(inst)

    path = clean.shortest_graphlike_error()
    print("effective distance:", len(path))
    for explained in circuit.explain_detector_error_model_errors(
        dem_filter=path, reduce_to_one_representative_error=True
    ):
        print(explained.dem_error_terms, "<-", explained.circuit_error_locations[0])
    return len(path)


if __name__ == "__main__":
    circuit = build_circuit()
    check_noiseless(circuit)
    print("noiseless check ok")
    effective_distance(circuit)
