# Surface

 ```Surface(dist: int, reps=1, errors=[0.0002,0.0041,.042,.012,.075]```

## Inputs
### reps : int
An integer value representing the desired number of repitions. Passes the number to my ```REPEAT {reps} {...``` block when initialising the circuit in STIM

### dist : int
An integer value representing the desired distance of the surface code. Note that dist must be an **odd number**, otherwise initialisatuion with raise `NotImplementedError`

### errors: list[float]
A list of floats representing the error rate for differernt components of the circuit. List is of the form

[single qubit, two qubit, measurement, idle, reset]

Note that errors passed as decimals not percentages. 

## Methods
### `Surface.write_to_stim()`
Returns the desired `stim.Circuit` object which represents the heavy hex surface code we want.

### `Surface.get_circuit()`
Returns `str` representation of surface with can be directly passed to `stim.Circuit()` to initilise cirucuit. This is used in the `Surface.write_to_stim()` but its nice to see the string for debugging or to put in your report idk.


# Notes
You should only need to import the `Surface` into your other projects, all other files are helper methods and classes for the `Surface`. Then you can use the methods above to initialise the stim circuit and pass it through your sims. 

I have found that the Scinter.Task initialisation takes some time - it may be becuase my solutions are kind of beefy for just initialising the circuit but see how you guys go.



