# Error details

Each channel is turned on or off, and given its own probability, in `error_model.json`. `ErrorModel` in `errors.py` reads that file. `CSSHoneycomb` loads it in `__init__` and appends the resulting Stim instructions while it builds the circuit string. A channel that is switched off, or that has no target qubits, adds nothing.

Stim draws the actual Pauli when the circuit is sampled. The circuit records the channel and its probability, not which Pauli fired in a given shot.

Each qubit gets one channel per layer. Qubits that took the layer's gate receive that gate's noise. Everyone else receives idle noise, from `_idle_keys` in `honeycomb.py`.

## Initialisation error

Initialisation error is one fault on a qubit that was just prepared, with probability `p_initialisation`. The Pauli depends on the prepared state, because only a Pauli that takes that state to an orthogonal one is a real fault.

It is implemented by `ErrorModel.initialisation_error`, using the `initialisation_error` flag. `get_measurement` calls it after the ancilla `R` and before `tick`. Those ancillas are always reset into `|0⟩`, so they take `X_ERROR`. `prepare_qubits` calls it after the data-qubit resets, and again after the Bell pair for the T qubit.

- `|0⟩`, from `R`: `X_ERROR`. This prepares `|1⟩` instead. A Y error is the same state.
- `|+⟩`, from `RX`: `Z_ERROR`. This prepares `|-⟩` instead. An X error does nothing to `|+⟩`.
- `|+i⟩`, from `RY`: `X_ERROR`. X and Z both take `|+i⟩` to `|-i⟩`. A Y error leaves `|+i⟩` unchanged.
- T data qubit: `DEPOLARIZE1` after the Bell pair, not after the `|0⟩` reset. Before the pair the qubit is still `|0⟩`, so X and Y are the same fault. After `H` on the T ancilla and `CNOT` onto the T data qubit, X, Y, and Z flip different correlations with the ancilla. Each Pauli has probability `p_initialisation / 3`. The Bell gates themselves are noiseless, and the T ancilla is left out of every channel.

## Hadamard error

Hadamard error is the chance that a qubit suffers a random X, Y, or Z after a Hadamard. Each Pauli has probability `p_hadamard / 3`.

It is implemented by `ErrorModel.hadamard_error`, which writes `DEPOLARIZE1(p_hadamard)` using the `hadamard_error` flag. `get_measurement` calls it after both Hadamard blocks in an X round: the block that puts ancillas into the X basis, and the block that rotates them back. A Z round has no Hadamard, so it has no Hadamard error.

## CNOT error

CNOT error is the chance that a CNOT applies a random two-qubit Pauli to its control and target. Each of the 15 non-identity Paulis has probability `p_cnot / 15`.

It is implemented by `ErrorModel.cnot_error`, which writes `DEPOLARIZE2(p_cnot)` using the `cnot_error` flag. The target list is the same control-target order as the `CNOT` instruction, so Stim pairs the same qubits. `get_measurement` calls it after each of the two CNOT rounds, before `tick`.

## Measurement error

Measurement error is the chance that the recorded measurement bit is flipped.

It is implemented by `ErrorModel.measurement_error`, which writes `X_ERROR(p_measurement)` before a Z-basis measurement and `Z_ERROR(p_measurement)` before an X-basis measurement, using the `measurement_error` flag. `get_measurement` calls it immediately before `M` and `MX` / `MZ`. Ancillas are always measured with `M`, including in an X round, because the Hadamard has already changed their basis, so they take `X_ERROR`. A boundary `MZ` takes `X_ERROR`. A boundary `MX` takes `Z_ERROR`.

## Idling error

Idling error is the possibility of a random X, Y, or Z on qubits that are not involved in any operation or measurement during a time step. Each Pauli has probability `p_idling / 3`.

It is implemented by `ErrorModel.idling_error`, which writes `DEPOLARIZE1(p_idling)` using the `idling_error` flag. The targets come from `CSSHoneycomb._idle_keys`, which returns every data qubit and ancilla that did not receive that layer's gate noise. `get_measurement` calls it in the reset layer, both Hadamard layers, both CNOT layers, and the measurement layer. `prepare_qubits` calls it on the ancillas, which are not reset in that step.

X rounds contain two Hadamard layers that Z rounds do not, so data qubits have two extra chances to idle in an X round. Z rounds are left shorter. No idle-only layers are added to pad them out.
