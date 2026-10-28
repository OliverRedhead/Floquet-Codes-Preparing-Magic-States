"""
Phenomenological-noise version of CSSHoneycomb with native two-body measurements.

Each round is a single layer:

    DEPOLARIZE1(p_data)  on every data qubit
    MPP(p_meas) ...      native XX / ZZ edge measurements (+ single-body boundary ones)
    TICK

There are no edge ancillas, no CNOTs, no Hadamards and no idling noise. The
only noise is:
    - p_data   : depolarizing noise on each data qubit, once per round, before measuring
    - p_meas   : classical flip probability on every MPP outcome
    - p_init   : preparation error on reset data qubits (X_ERROR for |0>, |+i>; Z_ERROR for |+>)
    - p_inject : DEPOLARIZE1 on the T data qubit after the Bell pair is made
The T ancilla (key 0) is noiseless everywhere.

Everything else (measurement bookkeeping, plaquette detectors, init detectors,
final readout, logical observables) is inherited unchanged.
"""

import stim

from honeycomb import CSSHoneycomb   # adjust to wherever CSSHoneycomb lives


class PhenomenologicalNoise:
    """
    Drop-in replacement for ErrorModel, exposing only what `prepare_qubits`
    needs (`initialisation_error`, `idling_error`) plus phenomenological helpers.
    A probability of 0 emits no noise instruction at all, so
    `PhenomenologicalNoise()` gives the noiseless circuit.
    """

    def __init__(self, p_data: float = 0.0, p_meas: float = 0.0,
                 p_init: float | None = None, p_inject: float | None = None,
                 *, T_ancilla_key: int = CSSHoneycomb.T_ANCILLA_KEY) -> None:
        self.p_data = p_data
        self.p_meas = p_meas
        self.p_init = p_data if p_init is None else p_init
        self.p_inject = self.p_init if p_inject is None else p_inject
        self.T_ancilla_key = T_ancilla_key

    def _channel(self, name: str, p: float, keys) -> str:
        keys = [k for k in keys if k != self.T_ancilla_key]
        if p <= 0 or not keys:
            return ""
        return f"{name}({p}) " + " ".join(map(str, keys)) + "\n"

    # --- interface used by CSSHoneycomb.prepare_qubits ---

    def initialisation_error(self, z_basis=(), x_basis=(), y_basis=(), t_basis=()) -> str:
        """|0> and |+i> get X_ERROR, |+> gets Z_ERROR, the T data qubit gets DEPOLARIZE1."""
        s = self._channel("X_ERROR", self.p_init, list(z_basis) + list(y_basis))
        s += self._channel("Z_ERROR", self.p_init, x_basis)
        s += self._channel("DEPOLARIZE1", self.p_inject, t_basis)
        return s

    def idling_error(self, *args, **kwargs) -> str:
        """No idling noise in the phenomenological model."""
        return ""

    # --- phenomenological helpers ---

    def data_error(self, keys) -> str:
        return self._channel("DEPOLARIZE1", self.p_data, keys)

    def mpp(self, targets: list[str]) -> str:
        """One MPP instruction with classical outcome-flip probability p_meas."""
        if not targets:
            return ""
        flip = f"({self.p_meas})" if self.p_meas > 0 else ""
        return f"MPP{flip} " + " ".join(targets) + "\n"


class PhenomenologicalHoneycomb(CSSHoneycomb):
    """CSSHoneycomb where every round is native two-body (and single-body) MPP."""

    def __init__(self, ncols: int = 5, nrows: int = 4,
                 noise: PhenomenologicalNoise | None = None, *,
                 config_filepath: str | None = None,
                 error_model: str | None = None) -> None:
        # The base constructor still loads an ErrorModel; we replace it right after.
        super().__init__(ncols, nrows, error_model, config_filepath=config_filepath)
        self.error_model = noise if noise is not None else PhenomenologicalNoise()
        self.error_model.T_ancilla_key = CSSHoneycomb.T_ANCILLA_KEY

    def initialise_circuit(self, coordinates='square') -> str:
        """Declare data qubits only; there are no edge ancillas in this model."""
        if coordinates != 'square':
            return super().initialise_circuit(coordinates)
        string = "# data qubits\n"
        for v in self.vertices:
            string += str(v.qubit) + "\n"
        return string

    def get_measurement(self, colour: str, flavour: str) -> str:
        """
        Measure every edge of `colour` as a native two-body `flavour flavour`
        Pauli, plus the single-body boundary measurements, in one MPP layer.

        Measurement order in the MPP instruction (edges, then boundary vertices)
        matches the order they are appended to `self.measurements`, so the
        rec[] lookups in `write_detectors` stay correct.
        """
        target_edges = [e for e in self.edges if e.colour == colour]
        # name-mangled private helpers of the base class
        target_vertices = self._CSSHoneycomb__get_boundary_targets(colour, flavour)

        string = f"# {colour}{flavour} round: native two-body measurements\n"

        data_keys = [v.key for v in self.vertices if v.key != CSSHoneycomb.T_ANCILLA_KEY]
        string += self.error_model.data_error(data_keys)

        targets = [
            f"{flavour}{e.vertices[0].key}*{flavour}{e.vertices[1].key}"
            for e in target_edges
        ]
        targets += [f"{flavour}{v.key}" for v in target_vertices]
        string += self.error_model.mpp(targets)

        edge_meas, vertex_meas = self._CSSHoneycomb__record_measurements(
            colour, flavour, target_edges, target_vertices
        )

        detectors = self.get_detectors()
        if self.isg_clock <= self._init_last_round:
            detectors += self.find_init_detectors()
        detectors += self._find_measurement_init_detectors(flavour, edge_meas, vertex_meas)
        string += self.write_detectors(detectors, self.clock)

        self._log_measured_ops(flavour, edge_meas, vertex_meas)
        self.isg_tick()

        # Tick at the end so every round gets a distinct measurement time
        # (`_readout_measurements` relies on that).
        string += self.tick()
        return string

    def _rounds(self, cycle) -> str:
        """Run get_measurement for each token ("rZ", "gX", ...) in order."""
        string = ""
        for token in cycle:
            colour = CSSHoneycomb._COLOUR_FROM_LETTER[token[0]]
            string += self.get_measurement(colour, token[-1])
        return string

    @classmethod
    def build_injection_circuit(
        cls,
        noise: PhenomenologicalNoise,
        *,
        ncols: int = 5,
        nrows: int = 3,
        config_filepath: str,
        ideal_cycle=("bZ", "rX", "gZ", "bX"),
        noisy_cycle=("rZ", "gX", "bZ", "rX", "gZ", "bX"),
    ) -> stim.Circuit:
        """
        Ideal injection, then one noisy cycle, then a noiseless readout.

        The ideal rounds prepare the lattice and inject the T state with no
        noise at all. The noisy cycle then runs with `noise.p_data` and
        `noise.p_meas` only; preparation noise (p_init, p_inject) stays off.
        Logical Z and X observables are added after the noisy cycle, followed
        by the noiseless final readout.
        """
        sim = cls(ncols, nrows, PhenomenologicalNoise(), config_filepath=config_filepath)
        sim.set_measurement_cycle([*ideal_cycle, *noisy_cycle])

        text = sim.initialise_circuit()
        text += sim.prepare_qubits()
        text += sim._rounds(ideal_cycle)

        sim.error_model = PhenomenologicalNoise(
            noise.p_data, noise.p_meas, p_init=0.0, p_inject=0.0
        )
        text += sim._rounds(noisy_cycle)

        last = noisy_cycle[-1]
        noiseless = str(stim.Circuit(text).without_noise())
        text += sim.write_logical_observable(
            noiseless, CSSHoneycomb._COLOUR_FROM_LETTER[last[0]], last[-1],
            bases=("Z", "X"),
        )
        text += sim.final_readout()
        return stim.Circuit(text)


def _build_rounds(sim: PhenomenologicalHoneycomb, cycle: list[str], initialisation) -> str:
    sim.set_measurement_cycle(cycle)
    text = sim.initialise_circuit()
    text += sim.prepare_qubits(initialisation)
    for colour, flavour in sim.measurement_cycle:
        text += sim.get_measurement(colour, flavour)
    return text


def build_phenomenological_circuit(
    cycle: list[str],
    *,
    noise: PhenomenologicalNoise,
    ncols: int = 5,
    nrows: int = 4,
    config_filepath: str | None = None,
    initialisation: dict | None = None,
    observables: bool = True,
    final_readout: bool = True,
    validate: bool = True,
) -> stim.Circuit:
    """
    Build the full noisy circuit for a measurement cycle such as
    ["rZ", "gX", "bZ", ...].

    A noiseless twin is built first. It is used to (a) solve for the
    OBSERVABLE_INCLUDE records, which needs a noiseless circuit, and (b) check
    that every detector and observable is deterministic. The observables are
    placed right after the last round, before the noiseless final readout.
    """
    def new(n):
        return PhenomenologicalHoneycomb(ncols, nrows, n, config_filepath=config_filepath)

    ideal = new(PhenomenologicalNoise())
    ideal_text = _build_rounds(ideal, cycle, initialisation)

    obs = ""
    if observables:
        colour, flavour = ideal.measurement_cycle[-1]
        obs = ideal.write_logical_observable(ideal_text, colour, flavour)
    tail = ideal.final_readout() if final_readout else ""

    if validate:
        # Raises if any detector / observable is non-deterministic.
        stim.Circuit(ideal_text + obs + tail).detector_error_model()

    noisy = new(noise)
    noisy_text = _build_rounds(noisy, cycle, initialisation)
    noisy_text += obs
    if final_readout:
        noisy_text += noisy.final_readout()

    return stim.Circuit(noisy_text)


if __name__ == "__main__":
    circuit = build_phenomenological_circuit(
        ["rZ", "gX", "bZ", "rX", "gZ", "bX"] * 2,
        noise=PhenomenologicalNoise(p_data=1e-3, p_meas=1e-3),
        config_filepath="config.json",
    )
    print(circuit.num_qubits, circuit.num_detectors, circuit.num_observables)
    print(circuit.detector_error_model(decompose_errors=False).num_errors)