from qubit import Qubit
import numpy as np

"""
Detectors in stim are set on measurements, stim holds measurements as just an ordered list
    where the last (index -1) item is the most recent measurement. This makes it difficult for 
    us to set up detectors, so a solution I found helpful last sem was to hold my own record of
    all measurements and just look up measurements in that list.

This may slow our curcuits down, but for now this is my best idea. If we have measurement objects,
    we may as well have other gates and operations like repeats to make everything more readable.
    These classes may also help us apply errors and things more easily.
"""


class Operation():
    """
    Operation base class. We will break our circuit into discrete time steps.
    """

    def __init__(self, time: int) -> None:
        self.time = time

class Gate(Operation):

    """
    Gate base class. We will consider measurements a gate as it's sort of the same thing as far 
    as stim is concerned.
    """

    # this will hold the prefix for stim
    OPERATION = ""

    def __init__(self, target : Qubit, time: int) -> None:
        """
        Initialise a Gate. 

        Parameters
        ----------
        target : Qubit
            the target qubit to apply the gate
        time : int
            Discrete time value to apply operation
        """
        self.target = target
        super().__init__(time)

    """
    these dunder methods are mainly to help with debugging, should make 
    degugging console more readable
    """

    def __repr__(self) -> str:
        """Return a string representation of the qubit (for debug)"""
        return f"{self.OPERATION}(target={self.target.key}, time={self.time})"
    
    def __str__(self) -> str:
        """Return a string representation of the qubit"""
        return f"{self.OPERATION} {self.target.key}\n"
    
    def __eq__(self, other) -> bool:
        """
        Define a notion of equality.

        Two gates are equal if they are of the same type (reset, measurement, cnot ...)
        and they act on the same target qubit at the same time.
        """
        if isinstance(other, type(self)) and ((self.target == other.target) and (self.time == other.time)):
            return True
        return False
    
    def __lt__(self, other) -> bool:
        """
        Define a notion of less than.

        This gate is less than another gate if it was applied to our circuit in an earlier time step.
        - Note that in this case, two gates do not need to be of the same type (reset, measurement, ...)
            we are only concerned with the time step they are applied in.
        """
        if not isinstance(other, Gate):
            return NotImplemented

        return self.time < other.time

    def __hash__(self) -> int:
        """
        Since we redefine equality, we override hash. In case we want to use hashable properties,
        we should explicitly redefine hash following our definition of equality.
        """
        return hash((type(self).__name__, self.target.key, self.time))
    



"""specific gate classes"""

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

    def __init__(self, target: Qubit, time: int, flavour: str = "Z") -> None:
        
        if flavour not in ['X', 'Y', 'Z']:
            raise ValueError(f"measurement.flavour must be 'X', 'Y' or 'Z', not {flavour}")
        self.flavour = flavour
        self.OPERATION += flavour

        super().__init__(target, time)

    def __eq__(self, other):
        """Need __eq__ to capture equality of measurement flavour too"""
        return super().__eq__(other) and (self.flavour == other.flavour)
    
    def __hash__(self) -> int:
        """
        Since we redefine equality, we override hash. In case we want to use hashable properties,
        we should explicitly redefine hash following our definition of equality.
        """
        return hash((type(self).__name__, self.target.key, self.time, self.flavour))
        
        

class Cnot(Gate):
    
    """
    A two-qubit gate with a target and control qubit. As this is the only two qubit gate we care
        about I will basically hardcode the different dunder methods.
    """

    OPERATION = "CNOT"

    def __init__(self, control: Qubit, target: Qubit, time: int) -> None:
        super().__init__(target, time)
        self.control = control

    def __repr__(self) -> str:
        """Return a string representation of the qubit (for debug)"""
        return f"{self.OPERATION}({self.control.key} -> {self.target.key}), t={self.time}"
    
    def __str__(self) -> str:
        """Return a string representation of the qubit"""
        return f"{self.OPERATION} {self.control.key} {self.target.key}\n"
    
    def __eq__(self, other) -> bool:
        """Need __eq__ to capture equality of the control qubits too"""
        return super().__eq__(other) and (self.control == other.control)
    
    def __hash__(self) -> int:
        """
        Since we redefine equality, we override hash. In case we want to use hashable properties,
        we should explicitly redefine hash following our definition of equality.
        """
        return hash((type(self).__name__, self.target.key, self.time, self.control))


"""Detectors are a special (classical) type of operation that compare measurement outcomes"""
    
class Detector(Operation):

    """
    In stim, for decoding purposes, detectors are given space-time coordinates. So we will need to specify
    this for some lovely visuals (detslice diagrams). This is basically why we break our circuit into discrete
    time steps. 
    
    Detectors make our life much more difficult.
    """

    def __init__(self, pos: tuple[int, int], time: int, rec: list[int]) -> None:
        """
        Initialise a Detectors instance.

        Parameters
        ----------
        pos : tuple[int, int]
            We need to specify a position. Note that this may not work in hexagonal coordinates.
        time : int
            We need to specify temporal location too.
        rec : list[int]
            This will be a list of negative integers which point stim to the correct measurements
            to compare. This is the difficult bit that requires us to keep a list of measurements
            in our circuits.
        """
        self.pos = pos
        self.rec = rec
        super().__init__(time)

    def __str__(self) -> str:
        """Let's us send a detector instance straight into stim."""
        x,y = self.pos
        string = f"DETECTOR({x},{y},{self.time})"
        for r in self.rec:
            string += f" rec[{r}]"
        return string