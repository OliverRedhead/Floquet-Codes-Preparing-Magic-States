import numpy as np
from qubit import Qubit


class Operation:

    def __init__(self, time: int) -> None:
        self.time = time


class Measurement(Operation):

    def __init__(self, time: int, target: Qubit, flavour: str) -> None:
        super().__init__(time)
        self.target = target
        self.flavour = flavour