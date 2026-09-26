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


# TODO may add colour attribute to measurement

class Measurement(Gate):

    OPERATION = "M"

    def __init__(self, target: Qubit, time: int, flavour: str) -> None:
        
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

    def __init__(
        self, 
        open_meas: list[Measurement] | None, 
        close_meas: list[Measurement] | None, 
        colour: str, 
        *, 
        flavour: str | None = None,
        open_time: int | None = None, 
        close_time: int | None = None
    ) -> None:
        """
        Initialise a Detectors instance.

        Parameters
        ----------
        open : list[Measurement] | None
            list of measurements that create syndrome at start of detector cell            
        close : list[Measurement] | None
            list of measurements that create syndrome at end of detector cell
        colour : str
            colour of detector cell
        flavour : str | None = None
            optional. flavour of detector cell, if not specified it is calculated.
        open_time : int | None = None
            optional. open time of detector cell, if not specified it is calculated.
        close_time : int | None = None
            optional. close time of detector cell, if not specified it is calculated.
        """
        
        open_meas = open_meas if open_meas is not None else []
        close_meas = close_meas if close_meas is not None else []
        
        self.open = open_meas
        self.close = close_meas
        self.colour = colour
        
        if open_time is not None:
            close_time = open_time + 4
        elif close_time is not None:
            open_time = close_time - 4
            
            if open_time < 0:
                raise NotImplementedError("You need to catch these detectors cells that open before our code starts!")
            
        else:
            if not (open_meas or close_meas):
                raise ValueError("Must specify at least one: open_meas, close_meas, open_time, close_time")
            
            open_time = open_meas[0].time if open_meas else close_meas[0].time - 4
            close_time = close_meas[0].time if close_meas else open_meas[0].time + 4
            
            if close_time - open_time != 4:
                raise NotImplementedError(f"Time on measurements does not work with our code.\
                                          opened at {open_time} and closed at {close_time} is not allowed.")
        
        # set open and close time attribute
        self.open_time = open_time
        self.close_time = close_time

        # set flavour attribute
        if flavour is None:
            # should have handled the case where both of these are empty already
            if open_meas: 
                self.flavour = open_meas[0].flavour
            else:
                self.flavour = close_meas[0].flavour
        else:
            self.flavour = flavour
        
        # make sure all measurements have same flavour.
        for m in open_meas + close_meas:
            if m.flavour != self.flavour:
                raise NotImplementedError("Measurement cannot have different flavour to Detector.")

        super().__init__(open_time)
        
    def is_valid(self):
        """
        returns bool if this is a valid detector. I originally set this up just to check the validity, 
        but may as well make sure everything is working properly while we are checking things.
        
        NOTE This may slow our generation down a bit so might be worth removing once we get things running.
        Also could just make this an attribute that gets updated everytime we add new measurement.
        """
        # 1. must have same number of open and close measurements.
        if (len(self.open) != len(self.close)) or (len(self.open) == 0):
            return False
        
        # 2. Measurements all have the same flavour.
        for m in self.open + self.close:
            if m.flavour != self.flavour:
                raise NotImplementedError("Measurement cannot have different flavour to Detector.")
    
        # 3. Measurements in open and close must occur at same respective non-negative time-steps
        for m in self.open:
            if m.time != self.open_time:
                raise NotImplementedError(f"Measurements in detector open must all be taken at the same time step.")
            if m.time < 0:
                return False
            
        for m in self.close:
            if m.time != self.close_time:
                raise NotImplementedError(f"Measurements in detector close must all be taken at the same time step.")
            if m.time < 0:
                raise ValueError("Somehow we have set up a measurement that closes before time starts.")
        
        
    def add_measurement(self, m: Measurement, end: str):
        """add a measurement to this detector

        Parameters
        ----------
        m : Measurement
            measurement object to add to this detector.
        end : str
            which end of the detector cell to add the measurement. 'open' or 'close'.
        """
        
        if m.flavour != self.flavour:
            raise NotImplementedError("Measurement cannot have different flavour to Detector.")
        
        if end == "open":
            if m.time != self.open_time:
                raise NotImplementedError(f"Measurements in detector open must all be taken at the same time step.")
            
            self.open.append(m)
            
        elif end == "close":
            if m.time != self.close_time:
                raise NotImplementedError(f"Measurements in detector close must all be taken at the same time step.")
            
            self.close.append(m)
            
        else:
            raise NotImplementedError(f"end must be 'open' or 'close', not {end}")
        
        
    def get_pos(self):
        """Compute the detector position as the COM of its qubits."""
        # TODO make sure this works
        return np.mean([m.target.pos for m in self.open], axis=0)
            