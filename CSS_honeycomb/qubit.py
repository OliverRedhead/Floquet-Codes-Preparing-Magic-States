class Qubit:
    """
    Qubit base class for my stim wrapper. We will use qubit instances in containers to set up our circuit.
    This base class means that we can hold the position and key of each qubit within itself and use some helper methods too.
    """

    def __init__(self, key: int, pos: tuple[float | int, float | int]) -> None:
        self.key = key        
        self.pos = pos

    
    def __repr__(self) -> str:
        """Return a string representation of the qubit (for debug)"""
        return f"Qubit(position={self.pos}, key={self.key})"
    
    def __iter__(self):
        """Iterate over the qubit's position and key."""
        x, y = self.pos
        yield from (self.key, x, y)

    def __str__(self) -> str:
        """Return a string representation of the qubit."""
        i,x,y = self
        return f"QUBIT_COORDS({x},{y}) {i}"
    
    def __eq__(self, other: object) -> bool:
        """
        Return True if this qubit is equal to another qubit.
        i.e. if this qubit is in same position with same key
        """
        if not isinstance(other, Qubit):
            return NotImplemented
        
        return (self.pos[0] == other.pos[0]
                and self.pos[1] == self.pos[1]
                and self.key == other.key)

    def __gt__(self, other:object) -> bool:
        """assert if this qubit is "greater than" another qubit - compare indices"""
        if not isinstance(other, Qubit):
            raise NotImplementedError
        if self.key > other.key:
            return True
        return False
    
    def __lt__(self, other:object) -> bool:
        """assert if this qubit is "less than" another qubit - compare indices"""
        if not isinstance(other, Qubit):
            raise NotImplementedError
        if self.key < other.key:
            return True
        return False
    
    def __hash__(self) -> int:
        """redefine hash after it is overriden by __eq__"""
        return hash((type(self).__name__, self.key, self.pos))

    