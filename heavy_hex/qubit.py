import numpy as np

class Qubit:

    def __init__(self, position: np.ndarray, index: int) -> None:
        self.position = position
        self.index = index
    
    def __eq__(self, other_qubit: object) -> bool:
        """
        Return True if this qubit is equal to another qubit.
        i.e. if this qubit is in same position with same index
        """
        if not isinstance(other_qubit, Qubit):
            return NotImplemented
        return (self.position[0] == other_qubit.position[0]
                and self.position[1] == self.position[1]
                and self.index == other_qubit.index)
    
    def __repr__(self) -> str:
        """Return a string representation of the qubit (for debug)"""
        return f"Qubit(position={self.position}, index={self.index})"
    
    def __iter__(self):
        """Iterate over the qubit's position and index."""
        yield from (self.position[0], self.position[1], self.index)

    def __str__(self) -> str:
        """Return a string representation of the qubit."""
        x,y,i = self
        return f"QUBIT_COORDS({x},{y}) {i}"

    def __gt__(self, other:object) -> bool:
        """assert if this qubit is "greater than" another qubit - compare indices"""
        if not isinstance(other, Qubit):
            raise NotImplementedError
        if self.index > other.index:
            return True
        return False
    
    def __lt__(self, other:object) -> bool:
        """assert if this qubit is "less than" another qubit - compare indices"""
        if not isinstance(other, Qubit):
            raise NotImplementedError
        if self.index < other.index:
            return True
        return False
    

class Data(Qubit):
    
    def __init__(self, position: np.ndarray, index: int) -> None:
        super().__init__(position, index)

    def __repr__(self) -> str:
        return f"Data(position={self.position}, index={self.index})"
    

class Syndrome(Qubit):

    def __init__(self, position: np.ndarray, index: int) -> None:
        super().__init__(position, index)
        self.measurements = []

    def __repr__(self) -> str:
        return f"Syndrome(position={self.position}, index={self.index})"
    
    def get_measurements(self):
        return self.measurements


class SyData(Data, Syndrome):

    def __init__(self, position: np.ndarray, index: int) -> None:
        super().__init__(position, index)
    
    def __repr__(self) -> str:
        return f"SyData(position={self.position}, index={self.index})"

class Bridge(Qubit):

    def __init__(self, position: np.ndarray, index: int) -> None:
        super().__init__(position, index)

    def __repr__(self) -> str:
        return f"Bridge(position={self.position}, index={self.index})"

