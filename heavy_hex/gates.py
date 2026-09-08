from rotated_bridge.qubit import Qubit, Data, Syndrome, Bridge, SyData
import numpy as np

class Operation():
    def __init__(self, time:int) -> None:
        self.time = time

class Gate(Operation):

    OPERATION = ""

    def __init__(self, target : Qubit, time: int) -> None:
        self.target = target
        self.time = time

    def __repr__(self) -> str:
        """Return a string representation of the qubit (for debug)"""
        return f"{self.OPERATION}(target={self.target.index}, time={self.time})"
    
    def __str__(self) -> str:
        """Return a string representation of the qubit"""
        return f"{self.OPERATION} {self.target.index}\n"
    
    def get_target_index(self) -> int:
        """getter method: returns index of target Qubit"""
        return self.target.index
    
    def get_target(self) -> Qubit:
        """getter method: returns target Qubit"""
        return self.target

    def isvalid(self) -> bool:
        if isinstance(self.target, Qubit):
            return True
        return False

    def __eq__(self, other) -> bool:
        if isinstance(other, Gate) and ((self.target == other.target) and (self.time == other.time)):
            return True
        return False


class Reset(Gate):

    OPERATION = "R"

    def __init__(self, target: Qubit, time: int) -> None:
        super().__init__(target, time)

class Hadamard(Gate):

    OPERATION = "H"

    def __init__(self, target: Qubit, time: int) -> None:
        super().__init__(target, time)

class Measurement(Gate):

    OPERATION = "M"
    count : int = -1

    def __init__(self, target: Qubit, time: int, xz: str) -> None:
        super().__init__(target, time)
        self.xz = xz
        if isinstance(target,Qubit):
            self.position = target.position
        if isinstance(target, Syndrome): 
            target.measurements.append(self)

    def set_count(self, count: int) -> None:
        self.count = count
    
    def get_count(self) -> int:
        return self.count
    
    def __iter__(self):
        yield from (self.position[0], self.position[1], self.time, self.count)

    def __repr__(self) -> str:
        """Return a string representation of the qubit (for debug)"""
        return f"{self.OPERATION}(target={self.target.index}, time={self.time}, xz={self.xz})"


class Cnot(Gate):

    OPERATION = "CNOT"

    def __init__(self, control: Qubit, target: Qubit, time: int) -> None:
        super().__init__(target, time)
        self.control = control

    def __repr__(self) -> str:
        """Return a string representation of the qubit (for debug)"""
        return f"{self.OPERATION}({self.control.index}->{self.target.index}), t={self.time}"
    
    def __str__(self) -> str:
        """Return a string representation of the qubit"""
        return f"{self.OPERATION} {self.control.index} {self.target.index}\n"
    
    def get_control_index(self) -> int:
        """getter method: returns index of control Qubit"""
        return self.control.index
    
    def get_indices(self) -> tuple[int, int]:
        """getter method: returns indices of control and target indices"""
        return self.control.index, self.target.index
    
    def get_control(self) -> Qubit:
        """getter method: returns control Qubit"""
        return self.control
    
    def get_qubits(self) -> tuple[Qubit, Qubit]:
        """getter method: returns control and target Qubit"""
        return self.control, self.target

    def __eq__(self, other) -> bool:
        if isinstance(other, Cnot) and ((super().__eq__(other)) and (self.control == other.control)):
                return True
        return False

    def isvalid(self) -> bool:
        if super().isvalid() and isinstance(self.control, Qubit):
            return True
        return False


class RepeatOpen(Operation):

    def __init__(self, n: int, time: int) -> None:
        super().__init__(time)
        self.n = n

    def open(self) -> str:
        string = "\nREPEAT " + str(self.n) + " {\n"
        return string
    
class RepeatClose(Operation):

    def __init__(self, time:int) -> None:
        super().__init__(time)
    
    def close(self) -> str:
        return "\n}\n"
    

class Detector(Operation):

    def __init__(self, position:np.ndarray, time: int, rec: list[int]) -> None:
        super().__init__(time)
        self.position = position
        self.rec = rec

    def __str__(self) -> str:
        x,y = self.position[0], self.position[1]
        string = f"DETECTOR({x},{y},{self.time})"
        for r in self.rec:
            string += f" rec[{r}]"
        return string