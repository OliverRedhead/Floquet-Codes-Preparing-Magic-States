import json
from pathlib import Path


class ErrorModel:
    """Stim noise channels loaded from error_model.json.

    Each channel has its own probability and an on/off flag. A disabled
    channel, or a channel with no target qubits, contributes an empty string.
    
    We specify the T_ancilla_key as a qubit to never add errors to
    """

    REQUIRED_KEYS = (
        "p_initialisation",
        "initialisation_error",
        "p_hadamard",
        "hadamard_error",
        "p_cnot",
        "cnot_error",
        "p_measurement",
        "measurement_error",
        "p_idling",
        "idling_error",
    )

    def __init__(self, config: dict, *, T_ancilla_key: int | None = None) -> None:
        missing = [key for key in self.REQUIRED_KEYS if key not in config]
        if missing:
            raise KeyError(f"error model config is missing keys: {missing}")
        self.config = config
        self.T_ancilla_key = T_ancilla_key

    @classmethod
    def load(cls, path: str | None = None, *, T_ancilla_key : int | None = None) -> "ErrorModel":
        """Load an error model. Defaults to error_model.json beside this file."""
        config_path = Path(path) if path is not None else Path(__file__).with_name("error_model.json")
        with config_path.open() as handle:
            return cls(json.load(handle), T_ancilla_key=T_ancilla_key)

    def _channel(self, instruction: str, probability_key: str, enabled_key: str, qubits) -> str:
        """Format one Stim noise instruction, or return an empty string."""
        if not self.config[enabled_key] or not qubits:
            return ""
        probability = self.config[probability_key]
        targets = " ".join(str(qubit) for qubit in qubits if qubit != self.T_ancilla_key)
        return f"{instruction}({probability}) {targets}\n"

    def initialisation_error(self, qubits) -> str:
        """X_ERROR after a reset into |0>."""
        return self._channel("X_ERROR", "p_initialisation", "initialisation_error", qubits)

    def hadamard_error(self, qubits) -> str:
        """DEPOLARIZE1 on qubits that just received a Hadamard."""
        return self._channel("DEPOLARIZE1", "p_hadamard", "hadamard_error", qubits)

    def cnot_error(self, pair_qubits) -> str:
        """DEPOLARIZE2 on CNOT pairs, in the same control-target order as the CNOT line."""
        return self._channel("DEPOLARIZE2", "p_cnot", "cnot_error", pair_qubits)

    def measurement_error(self, z_basis, x_basis) -> str:
        """Bit-flip noise before a measurement: X before Z-basis, Z before X-basis."""
        return (
            self._channel("X_ERROR", "p_measurement", "measurement_error", z_basis)
            + self._channel("Z_ERROR", "p_measurement", "measurement_error", x_basis)
        )

    def idling_error(self, qubits) -> str:
        """DEPOLARIZE1 on qubits that took no gate in this layer."""
        return self._channel("DEPOLARIZE1", "p_idling", "idling_error", qubits)
