import json

import matplotlib.pyplot as plt
import numpy as np
import stim

from qubit_containers import Vertex, Edge, Plaquette, BoundaryPlaquette
from operations import Measurement, PauliString


class CSSHoneycomb:
    """
    Owns the vertices, edges and plaquettes of a lattice, keyed by key.
    Responsible for construction, lookup, incidence queries and validation.
    """

    """
    initialisation methods
    """
    
    T_ANCILLA_KEY = 0

    def __init__(self, ncols : int = 5, nrows : int = 4, *, config_filepath: str | None = None) -> None:
        """
        Initialise CSSHoneycomb instance. Automatically generates vertices (and qubits) based on input nrows and ncols.

        Parameters
        ----------
        nrows : int = 4
            number of rows of qubits
        ncols : int = 5
            number of columns of qubits
        config_filepath: str | None
            filepath to config, default None.

        Notes
        -----
        - Coordinates are initialised in integer (square) system. This is more natural for the heavy hex. 
            I will try and make sure there are adequate methods to let us set up the circuit in both square 
            and hexagonal coordinates, but square makes most sense when initialising qubits and things.
        """

        self.nrows = nrows
        self.ncols = ncols

        self.vertices = self.__initialise_vertices()
        self.edges = self.__initialise_edges()
        self.plaquettes = self.__initialise_plaquettes()
        self.plaquettes += self.__initialise_boundary_plaquettes()
        
        self.qubit_coords = {q.pos : q.key for q in self.vertices}
        
        self.vertices.append(Vertex(0, (self.ncols//2, -1))) # ancilla for T state
        
        self.measurements = []
        
        # This counts our circuit level ticks
        self.clock = 0
        
        # This counts ISG updates
        self.isg_clock = 0
        
        # set config
        self.preparation = None
        self.logical_schedule = None   # {"rXX": {"X": [keys], "Z": [keys]}, ...}

        if config_filepath is not None:
            with open(config_filepath) as f:
                raw = json.load(f)

            self.preparation = {
                self.qubit_coords[tuple(q["pos"])]: q["init"] for q in raw["qubits"]
            }
            self.logical_schedule = {
                name: {
                    basis: [self.qubit_coords[tuple(pos)] for pos in ops[f"{basis}_logical"]]
                    for basis in ("X", "Z")
                }
                for name, ops in raw["logical_schedule"].items()
            }
        
    
    
    def __initialise_vertices(self) -> list[Vertex]:
        """Initialise vertices, note that data qubits are initialised in this step too"""

        vertices = [] 
        key = 1 # we will put our ancilla at key=0

        for y in range(self.nrows):    
            n_vertices = self.ncols if y in (0, self.nrows - 1) else self.ncols + 1

            for x in range(n_vertices):
                
                x_pos = x + y - 1 if y > 0 else x

                vertices.append(Vertex(key, (x_pos, y)))
                key += 1

        return vertices

    def __initialise_edges(self) -> list[Edge]:
        """
        Initialise edges.

        Must be called after vertices have been initialised.
        Each edge is defined by the two vertices it joins.
        """

        edges = []

        def make_edge(v0: Vertex, v1: Vertex, direction: str) -> Edge:
            edge_key = len(edges)
            colour = CSSHoneycomb.__colour_edge(v0.pos[0], direction)
            ancilla_key = len(self.vertices) + edge_key + 1 # +1 is because I add the ancilla for logicals later.

            return Edge(
                edge_key,
                v0=v0,
                v1=v1,
                colour=colour,
                ancilla_key=ancilla_key
            )

        for i, v0 in enumerate(self.vertices):
            x, y = v0.pos

            # Vertical edge
            if y < self.nrows - 1 and (x + y) % 2 == 0:
                v1 = self.vertices[i + self.ncols]
                edges.append(make_edge(v0, v1, "vertical"))

            # Horizontal edge
            if x < self.ncols + y - 1:
                # Last vertex of the bottom row has no horizontal neighbour
                if y == self.nrows - 1 and x == self.ncols + self.nrows - 3:
                    continue

                v1 = self.vertices[i + 1]
                edges.append(make_edge(v0, v1, "horizontal"))

        return edges

    def __initialise_plaquettes(self) -> list[Plaquette]:
        """
        Initialise plaquettes.

        Must be called after vertices have been initialised.
        Each plaquette is defined by the six vertices around it.
        """

        vertex_coords = {
            v.get_coordinates(): v
            for v in self.vertices
        }

        offsets = [
            (0, 0),
            (1, 0),
            (2, 0),
            (2, 1),
            (1, 1),
            (0, 1)
        ]

        plaquettes = []

        for (x, y) in vertex_coords:

            # Only consider one sublattice of possible plaquette origins
            if (x + y) % 2 == 1:
                continue

            positions = [
                (x + dx, y + dy)
                for dx, dy in offsets
            ]

            # Skip positions where the complete plaquette does not exist
            if not all(pos in vertex_coords for pos in positions):
                continue

            vertices = [
                vertex_coords[pos]
                for pos in positions
            ]

            colour = CSSHoneycomb.__colour_plaquette(x)

            plaquettes.append(
                Plaquette(
                    key=len(plaquettes),
                    colour=colour,
                    vertices=vertices
                )
            )

        return plaquettes

    def __initialise_boundary_plaquettes(self) -> list[BoundaryPlaquette]:
        """Initialise the plaquettes along all four boundaries.

        The top/bottom boundaries carry Z-flavoured plaquettes, the
        left/right boundaries carry X-flavoured ones. Each boundary is
        just a list of (x, y) anchor points, so they all share the same
        plaquette-building loop below.
        """
        vertex_coords = {
        vertex.get_coordinates(): vertex
        for vertex in self.vertices
        }

        offsets = [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1)]

        boundaries = [
            ("Z", [(x, -1) for x in range(1, self.ncols - 2, 2)]),                       # top
            ("Z", [(x, self.nrows - 1) for x in range(self.nrows - 1, self.ncols + self.nrows - 4, 2)]),  # bottom
            ("X", [(y - 2, y) for y in range(1, self.nrows - 1)]),                       # left
            ("X", [(self.ncols - 1 + y, y) for y in range(0, self.nrows - 2)]),               # right
        ]

        boundary_plaquettes = []

        for flavour, anchors in boundaries:
            for x, y in anchors:
                plaquette = CSSHoneycomb.__make_boundary_plaquette(
                    x, y, vertex_coords, offsets,
                    index=len(self.plaquettes) + len(boundary_plaquettes),
                    flavour=flavour,
                )
                if plaquette is not None:
                    boundary_plaquettes.append(plaquette)

        return boundary_plaquettes
        
    
    """
    stim interaction methods
    """

    # circuit initialisation

    def initialise_circuit(self, coordinates='square') -> str:
        """
        Converts our CSSHoneycomb into stim language.

        NOTE For now, the circuit only works in square coordinates.
        """
        
        if coordinates != 'square': # NOTE this doesnt work anymore.
            return self.__initialise_circuit_hex()
        
        string = "# data qubits\n"
        for v in self.vertices:
            string += str(v.qubit) + "\n"

        string += "# ancilla qubits\n"
        for e in self.edges:
            string += str(e.ancilla) + "\n"
        
        return string
    
    def __initialise_circuit_hex(self):
        string = ""        
        
        coords = [v.pos for v in self.vertices]
        
        hex_coords = CSSHoneycomb.square_to_hex(coords)
        for i, (x, y) in enumerate(hex_coords):
            string += f"QUBIT_COORDS({x}, {y}) {i}\n"
        return string

    def tick(self) -> str:
        """increase the clock and return the stim TICK command"""
        self.clock += 1
        return "\nTICK\n"
    
    def isg_tick(self) -> None:
        """increment isg clock"""
        self.isg_clock += 1
    
    # measurement protocol

    @staticmethod
    def __get_CNOTs(edge: Edge, flavour="Z", r=0) -> tuple[int, int]:
        """
        Helper method to get the CNOT schedule of a particular edge.

        Parameters
        ----------
        edge : Edge
            The edge we want to measure
        flavour : str = "Z"
            The flavour of measurement. This will determine the direction of CNOTs
        r : int = 0 | 1
            The 'round' of CNOTs we are applying. Should be binary. Determines if we
            apply CNOT to edge.vertex0 or edge.vertex1

        Returns
        -------
        cnot : tuple[int, int]
            A tuple of keys of the cnot. Should be in order (control, target).

        Notes
        -----
        Need to make sure that the order of this is correct!
        """
        vertex = edge.vertices[r]
        ancilla = edge.ancilla

        if flavour == "X":
            # ZZ: ancilla -> data
            return ancilla.key, vertex.key

        elif flavour == "Z":
            # XX: data -> ancilla
            return vertex.key, ancilla.key

        else:
            raise ValueError("flavour must be 'X' or 'Z'")

    def __get_boundary_targets(self, colour, flavour):
        """Helper to get all the single-body measurement targets on the boundary this round.

        Parameters
        ----------
        colour : str
            either 'red', 'green' or 'blue'. Determines boundary conditions.
        flavour : str
            Either 'X' or 'Z'. Determines boundary conditions.
        """
        
        vertex_coords = {
            v.pos : v
            for v in self.vertices
        }
        
        target_vertices = [
            v
            for p in self.plaquettes
            if isinstance(p, BoundaryPlaquette)
            and p.flavour == flavour
            and p.colour != colour
            for v in p.vertices.cells
        ] + [
            vertex_coords[(0,0)], 
            vertex_coords[(self.ncols + self.nrows - 3, self.nrows - 1)]
        ]
        
        targets = []
        for v in target_vertices:
            coboundary_colours = [e.colour for e in v.coboundary()]
            if colour not in coboundary_colours:
                targets.append(v)
            
        return list(set(targets))
    
    def _get_logical_keys(self, colour: str, flavour: str) -> dict[str, list[int]]:
        """
        Look up the data-qubit keys supporting each logical operator.

        The logical representative changes as the ISG evolves, so it is keyed by
        the round that has just been measured.

        Parameters
        ----------
        colour : str
            'red', 'green' or 'blue'. Colour of edges measured this round.
        flavour : str
            'X' or 'Z'. Flavour measured this round.

        Returns
        -------
        dict[str, list[int]]
            {"X": [...], "Z": [...]}: qubit keys of the X and Z logicals.
            The T ancilla is not included; see `_logical_pauli`.
        """
        if self.logical_schedule is None:
            raise ValueError("No logical_schedule loaded; pass config_filepath.")

        name = f"{colour[0]}{flavour}{flavour}"   # e.g. ('red', 'X') -> 'rXX'
        try:
            return self.logical_schedule[name]
        except KeyError:
            raise ValueError(f"No logical entry for round {name!r} in the schedule.") from None

    def __logical_pauli(self, colour: str, flavour: str, basis: str) -> dict[int, str]:
        """
        Build the Pauli support {qubit_key: 'X' | 'Z'} of a logical.

        The T ancilla is included because the logical is entangled with it
        through the Bell pair, so only the product with the ancilla's Pauli is
        fixed by the initial state.
        """
        if basis not in ("X", "Z"):
            raise ValueError(f"basis must be 'X' or 'Z', not {basis!r}")

        keys = self._get_logical_keys(colour, flavour)[basis]
        return {k: basis for k in (CSSHoneycomb.T_ANCILLA_KEY, *keys)}

    def __write_observable(self, circuit_so_far: str, index: int, pauli_by_key: dict[int, str]) -> str:
        """
        Emit an OBSERVABLE_INCLUDE for a Pauli operator, adding the measurement
        records that make its value deterministic in a dynamic code.

        In a dynamic code the logical representative changes between rounds, so
        the Pauli operator alone is not a deterministic observable: its value
        also depends on the outcomes of the checks multiplied in along the way.
        Stim's flow solver finds the measurement records whose parity must be
        included with the Pauli for the observable to be deterministic.

        Parameters
        ----------
        circuit_so_far : str
            Noiseless stim circuit text up to the point where the operator is
            evaluated. Record offsets are relative to the end of this circuit,
            so the returned line must be appended directly after it, before any
            further measurements are added.
        index : int
            Observable index, as in OBSERVABLE_INCLUDE(index).
        pauli_by_key : dict[int, str]
            Support of the operator as {qubit_key: 'X' | 'Z'}, e.g.
            {0: "Z", 2: "Z", 7: "Z"}.

        Returns
        -------
        str
            A single line of the form
            `OBSERVABLE_INCLUDE(index) Z0 Z2 Z7 rec[-5] rec[-12]\n`.

        Raises
        ------
        ValueError
            If no combination of measurements makes the operator a valid flow,
            meaning it is not a logical of the state at this point. This usually
            indicates a wrong `logical_schedule` entry for the last measured
            round, or a circuit that was not built up to that round.

        Notes
        -----
        - Pauli targets in OBSERVABLE_INCLUDE read the operator non-destructively,
        so nothing is measured and `self.measurements` is unaffected.
        - Run this on the noiseless circuit. The records found remain valid once
        noise instructions are inserted.
        - Requires stim >= 1.15 for `Circuit.solve_flow_measurements`.
        """
        circuit = stim.Circuit(circuit_so_far)

        logical = stim.PauliString(circuit.num_qubits)
        for k, p in pauli_by_key.items(): 
            logical[k] = p

        # Ask stim which measurements make "nothing -> logical" a valid flow.
        (recs,) = circuit.solve_flow_measurements([stim.Flow(output=logical)])
        if recs is None:
            raise ValueError(
                "Pauli string is not a logical of the state at this point. "
                "Check the logical_schedule entry for the last measured round."
            )

        n = circuit.num_measurements
        recs = [r - n if r >= 0 else r for r in recs]   # absolute -> rec[-k] offsets

        paulis = " ".join(f"{p}{k}" for k, p in pauli_by_key.items())
        rec_str = " ".join(f"rec[{r}]" for r in recs)
        return f"OBSERVABLE_INCLUDE({index}) {paulis} {rec_str}\n"

    def write_logical_observable(self, circuit_so_far: str, colour: str, flavour: str,
                                bases: tuple[str, ...] = ("Z", "X")) -> str:
        """
        Emit one OBSERVABLE_INCLUDE per requested logical basis, for the logical
        representatives valid after round (colour, flavour).

        Both logicals can be tracked in the same run because each is paired with
        the T ancilla's Bell partner Pauli (X_L X_anc and Z_L Z_anc commute).

        Parameters
        ----------
        circuit_so_far : str
            Noiseless stim circuit text up to and including the round. Must end
            right after that round's `get_measurement` output.
        colour, flavour : str
            The round that was just measured; selects the logical representatives.
        bases : tuple[str, ...] = ("Z", "X")
            Which logicals to emit. Observable index follows the order given,
            so the default gives Z -> OBSERVABLE 0 and X -> OBSERVABLE 1.

        Returns
        -------
        str
            One OBSERVABLE_INCLUDE line per basis.

        Notes
        -----
        All lines are computed against the same `circuit_so_far`. OBSERVABLE_INCLUDE
        adds no measurement records, so the rec[] offsets are the same for each.
        """
        return "".join(
            self.__write_observable(circuit_so_far, i, self.__logical_pauli(colour, flavour, b))
            for i, b in enumerate(bases)
        )
    
    """
    This dictionary holds the detector cell schedule. It is of the form
           dict[ tuple[str] : dict[ str : str ] ]
    The first dimension is to call on the measurement flavour and colour.
    For example, if we are measuring gXX this round, we might call
        syndrome = SYNDOME_SCHEDULE[('green', 'X')]
    this would give us another dictionary of the form
            dict[str : str].
    Now we can find what colour detector cells **open** and **close** with those keywords.
    """
    SYNDOME_SCHEDULE = {
        ('red', 'Z')    : {'open' : 'blue',  'close' : 'green'},
        ('green', 'X')  : {'open' : 'red',   'close' : 'blue'},
        ('blue', 'Z')   : {'open' : 'green', 'close' : 'red'},
        ('red', 'X')    : {'open' : 'blue',  'close' : 'green'},
        ('green', 'Z')  : {'open' : 'red',   'close' : 'blue'},
        ('blue', 'X')   : {'open' : 'green', 'close' : 'red'},
    }

    def __record_measurements(self, colour, flavour, target_edges, target_vertices):
        """Record measurements for bulk and boundary plaquettes.

        Parameters
        ----------
        colour : str
            'red', 'green' or 'blue'. Colour of edge measurements this round.
        flavour : str
            'X' or 'Z'. Flavour of edge measureremenets this round.
        target_edges : list[Edge]
            list of target edges measured this round, produced by get_measurement() method
        target_vertices : list[Vertex]
            list of target vertices measured this round, produced by get_measurement() method
        """

        open_colour = CSSHoneycomb.SYNDOME_SCHEDULE[(colour, flavour)]['open']
        close_colour = CSSHoneycomb.SYNDOME_SCHEDULE[(colour, flavour)]['close']

        # construct plaquettes in bulk with edges.
        for edge in target_edges:
            m = Measurement(edge.ancilla, self.clock, flavour)
            self.measurements.append(m)          # once per real M instruction

            for p in edge.coboundary():
                if p.colour == open_colour:
                    p.add_measurement(m, self.isg_clock, "open")
                elif p.colour == close_colour:
                    p.add_measurement(m, self.isg_clock, "close")

        # construct boundary plaquettes with single qubit measurements where applicable.

        for vertex in target_vertices:
            m = Measurement(vertex, self.clock, flavour)
            self.measurements.append(m)

            for edge in vertex.coboundary():
                for p in edge.coboundary():
                    if not isinstance(p, BoundaryPlaquette):
                        continue
                    if p.colour == open_colour:
                        p.add_measurement(m, self.isg_clock, "open")
                    elif p.colour == close_colour:
                        p.add_measurement(m, self.isg_clock, "close")
                                     
    def get_detectors(self):
        """Get a list of the valid detectors that close at the given time step."""
        detectors = []
        for p in self.plaquettes:
                    
            d = p.detector_closes(self.isg_clock)
            if d is not None:
                detectors.append(d)
            
        return detectors
    
    def write_detectors(self, detectors, time) -> str:
        """Given a list of valid detectors, write them in stim.
        
        Parameters
        ----------
        detectors : list[Detector]
            list of valid detectors.
            
        Returns 
        -------
        detector_string : str
            Detector string to pass to stim.
        """
        
        def find_measurement_index(m : Measurement):
            """find the rec index of a target measurement in the self.measurements list"""
            idx = self.measurements.index(m)
            return idx - len(self.measurements)
        
        detector_string = "\n"
        for d in detectors:
            x,y = d.get_pos()
            detector_string += f"DETECTOR({x}, {y}, {time}) "
            for m in d.open_meas + d.close_meas:
                idx = find_measurement_index(m)
                detector_string += f"rec[{idx}] "
            
            detector_string += "\n"
        
        return detector_string        

    def get_measurement(self, colour: str, flavour: str) -> str:
        """
        By specifying a colour and flavour of measurement, we can 
        measure all edges of the given colour in the basis specified.
        
        I want to include boundaries in here too.

        Parameters
        ----------
        colour : str
            either 'red', 'green' or 'blue'. Indicates what colour of edges we should measure.
        flavour : str
            Either 'X' or 'Z', indicates the measurement basis. Note that the 'X' circuit is a 
            little bit longer.

        Returns
        -------
        circuit string : str
            A string in `stim` syntax that we can pass directly in. Should be commented adequately too.
            
        Notes
        -----
        Our single qubits spend some time idling here. I'm not sure how avoidable this is.
        """
        string = ""

        target_edges = [e for e in self.edges if e.colour == colour]
        ancilla_keys = [e.ancilla.key for e in target_edges]
        
        target_vertices = self.__get_boundary_targets(colour, flavour)
        vertex_keys = [v.key for v in target_vertices]

        # Initialise ancillas in |0>
        string += f"# prepare ancillas in {flavour} basis\nR "
        for a in ancilla_keys:
            string += str(a) + " "
        string += self.tick()
        
        # If measuring X stabilizer, put ancillas in X basis
        if flavour == "X":
            string += "H "      # Hadamard
            for a in ancilla_keys:
                string += str(a) + " "                
            string += self.tick()

        # First round of CNOTs between v0 and ancilla
        string += "# fold stabilizers\nCNOT "
        for edge in target_edges:
            c, t = CSSHoneycomb.__get_CNOTs(edge, flavour=flavour, r=0)
            string += f"{c} {t} "
        string += self.tick()

        # Second round of CNOTs between v1 and ancilla
        string += "CNOT "
        for edge in target_edges:
            c, t = CSSHoneycomb.__get_CNOTs(edge, flavour=flavour, r=1)
            string += f"{c} {t} "
        string += self.tick()

        # If measuring X stabilizer, put ancillas back in Z basis
        if flavour == "X":
            string += "# rotate ancillas back to Z basis\nH "
            for k in ancilla_keys + vertex_keys:
                string += str(k) + " "
            string += self.tick()

        # measure ancillas
        string += f"#{colour}{flavour} edges\nM "
        for a in ancilla_keys:
            string += str(a) + " "
        
        # measure boundary data qubits (always in the Z basis)
        if vertex_keys:
            string += f"\n#{colour}{flavour} boundary\nM "
            for v in vertex_keys:
                string += str(v) + " "

            # undo the basis change so the qubits are back in the X basis
            if flavour == "X":
                string += "\n# restore boundary qubits to X basis\nH "
                for v in vertex_keys:
                    string += str(v) + " "
            
                
        self.__record_measurements(colour, flavour, target_edges, target_vertices) # NOTE not sure if I need colour here
        
        detectors = self.get_detectors()
        string += self.write_detectors(detectors, self.clock)
        
        self.isg_tick()
    
        return string

    def final_readout(self) -> str:
        """
        Perform a final, noiseless readout of every stabilizer still in the ISG.

        Finds every detector that has been opened but not yet closed, measures
        the corresponding plaquette stabilizer directly with a single MPP
        instruction, and closes the detector with that measurement. The DETECTOR
        lines for all of these detectors are then written, so that each one
        compares its opening edge measurements against the final readout.

        This should be called once, after the last round of `get_measurement`.
        Because a closed detector is skipped by `get_open_detectors`, a second
        call returns an empty string rather than re-measuring anything.

        No error model should be applied to this round.

        Returns
        -------
        circuit string : str
            A string in `stim` syntax: a TICK, a comment header, one
            `MPP` line per open detector (e.g. `MPP X3*X4*X9`), then the
            DETECTOR lines closing them. Returns an empty string if there are
            no open detectors.

        Notes
        -----
        - Each `PauliString` is appended to `self.measurements` before
        `write_detectors` is called, since the `rec[-k]` offsets are looked up
        from that list. The order of `MPP` lines must therefore match the order
        of appends, which it does because both happen in the same loop.
        - The detector is closed directly with `Detector.add_measurement`, not via
        `Plaquette.add_measurement`, so that the readout is attached to the
        specific open detector and not to an older, already-closed one of the
        same flavour on that plaquette.
        - Advances `self.clock` (via `tick`) and `self.isg_clock` (via
        `isg_tick`), like `get_measurement`.
        """
        open_pairs = [
            (p, d)
            for p in self.plaquettes
            for d in p.get_open_detectors()
        ]

        if not open_pairs:
            return ""

        string = self.tick()
        string += "# final readout: measure out remaining ISG stabilizers\n"

        for p, d in open_pairs:
            qubits = [v.qubit for v in p.vertices]   # same iteration as get_MPPstring

            m = PauliString(qubits, self.clock, d.flavour)
            self.measurements.append(m)
            d.add_measurement(m, "close")            # close this exact detector

            string += p.get_MPPstring(d.flavour) + "\n"

        string += self.write_detectors([d for _, d in open_pairs], self.clock)
        self.isg_tick()

        return string

    def prepare_qubits(self) -> str:
        """
        Prepare all qubits for the start of the circuit.

        Every qubit is reset in the Z basis by default. If a preparation config
        was loaded, qubits listed there are reset in the requested basis instead:
        'X' -> |+>, 'Y' -> |+i>, 'Z' (or '0') -> |0>.

        The qubit marked 'T' is reset to |0> and then entangled with the T
        ancilla (key 0) in a Bell pair: H on the ancilla, then CNOT
        ancilla -> T qubit. This gives (|00> + |11>)/sqrt(2).

        Returns
        -------
        circuit string : str
            Stim reset instructions, followed by the Bell pair circuit if a T
            qubit is configured, each layer ending in a TICK.
        """
        RESET_FOR_STATE = {"Z": "R", "X": "RX", "Y": "RY", "T": "R"}

        resets = {"R": [], "RX": [], "RY": []}
        t_key = None

        if self.preparation is None:
            resets["R"] = [v.key for v in self.vertices]
        else:

            t_keys = [k for k, s in self.preparation.items() if s == "T"]
            if len(t_keys) != 1:
                raise ValueError(f"Expected exactly one T qubit, found {len(t_keys)}")
            t_key = t_keys[0]

            for v in self.vertices:
                state = self.preparation.get(v.key, "Z")
                if state not in RESET_FOR_STATE:
                    raise ValueError(f"Unknown preparation state {state!r} on qubit {v.key}")
                resets[RESET_FOR_STATE[state]].append(v.key)

        string = "# initialise qubits\n"
        for instruction, keys in resets.items():
            if keys:
                string += f"{instruction} " + " ".join(map(str, keys)) + "\n"
        string += self.tick()

        # Bell pair between the T ancilla and the T qubit
        if t_key is not None:
            string += "# Bell pair: T ancilla <-> T qubit\n"
            string += f"H {CSSHoneycomb.T_ANCILLA_KEY}\n"
            string += f"CNOT {CSSHoneycomb.T_ANCILLA_KEY} {t_key}" + self.tick()

        return string
        

    """
    helper methods
    """

    @staticmethod
    def __colour_plaquette(x):
        """
        Helper for colouring plaquettes.
        """
        brg = ['blue', 'red', 'green']
        return brg[x%3]
    
    @staticmethod
    def __colour_edge(x, direction):
        """
        Helper for colouring edges. 
        """
        rgb = ['red', 'green', 'blue']
        if direction == "vertical":
            return rgb[(x + 1) % 3]
        return rgb[x % 3]

    @staticmethod
    def __make_boundary_plaquette(x, y, vertex_coords, offsets, index, flavour):
        """Build a single boundary plaquette anchored at (x, y).

        Returns None if none of the plaquette's vertices exist in
        `vertex_coords` (i.e. it would be entirely off-grid).
        """
        positions = [(x + dx, y + dy) for dx, dy in offsets]
        candidates = [vertex_coords.get(pos) for pos in positions]

        boundary_vertices = [v for v in candidates if v is not None]
        if not boundary_vertices:
            return None

        missing_positions = [
            pos for pos, v in zip(positions, candidates) if v is None
        ]

        return BoundaryPlaquette(
            index,
            colour=CSSHoneycomb.__colour_plaquette(x),
            vertices=boundary_vertices,
            missing_positions=missing_positions,
            flavour=flavour,
        )

    @staticmethod
    def square_to_hex(coords, scale=1):

        coords = np.asarray(coords, dtype=float)

        # Make sure we have an (N, 2) array
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError(
                f"coordinates must have shape (N, 2), got {coords.shape}"
            )

        a = np.sqrt(3) / 4

        x = coords[:, 0].copy()
        y = coords[:, 1].copy()

        even = (x % 2 == 0)

        y[even] += 2 * a * (y[even] // 2)

        y[~even] += (
            a * ((y[~even] + 1) // 2)
            + a * ((y[~even] - 1) // 2)
        )

        x *= 2 * a

        return np.column_stack((x, y)) * scale

    @staticmethod
    def hex_to_square(coords, scale=1):
        """
        Inverse of square_to_hex: given hexagonal coordinates, recover the
        original integer (x, y) square-grid coordinates.
        """
        square_coords = []
        a = np.sqrt(3) / 4

        for xh, yh in coords:
            xh, yh = xh / scale, yh / scale
            y = int(round(yh / (2 * a)))

            # continuous estimate of x, then verify/correct against the
            # exact forward formula (handles the floor-division steps)
            x_est = int(round(xh / (1 + a)))
            x = None
            for cand in range(x_est - 2, x_est + 3):
                if y % 2 == 0:
                    test = cand + (2*a) * (cand // 2)
                else:
                    test = cand + a * ((cand + 1) // 2) + a * ((cand - 1) // 2)
                if np.isclose(test, xh):
                    x = cand
                    break

            if x is None:
                raise ValueError(f"Could not invert hex coordinate ({xh}, {yh})")

            square_coords.append((x, y))

        return np.array(square_coords)
    

    """
    visualisation methods
    """

    PALETTE = {
        'red' : '#e74c3c', 
        'green' : '#2ecc71',
        'blue' : '#3498db', 
    }

    def plot_surface(self, coordinates="square", bulge=0.35, n_arc=16, savepath = None, show=True):

        """Visualise the surface."""

        if coordinates == "square":

            def coord_map(x):
                return np.asarray(x, dtype=float)

        elif coordinates == "hex":

            coord_map = self.square_to_hex  # type: ignore

        else:

            raise ValueError(
                f"`coordinates` must be 'square' or 'hex', not {coordinates!r}"
            )

        for p in self.plaquettes:

            p_coords = coord_map(
                np.asarray(
                    p.get_coordinates(ordered=True),
                    dtype=float
                )
            )

            coords = np.array([v.get_coordinates() for v in self.vertices])
            mapped_coords = coord_map(
               np.asarray(coords, dtype=float)
            )

            if len(p) == 6:

                plt.fill(
                    p_coords[:, 0],
                    p_coords[:, 1],
                    alpha=0.3,
                    facecolor=p.colour,
                    edgecolor="black",
                    zorder=1
                )

            elif len(p) in [3, 4]:

                gap_index = p.get_boundary_gap_index()

                if gap_index is None or not p.missing_positions:
                    patch_coords = p_coords
                else:
                    missing_mapped = coord_map(
                        np.asarray(p.missing_positions, dtype=float)
                    )
                    missing_centre = missing_mapped.mean(axis=0)
                    p_centre = p_coords.mean(axis=0)

                    outward_dir = missing_centre - p_centre
                    norm = np.linalg.norm(outward_dir)
                    outward_dir = (
                        outward_dir / norm if not np.isclose(norm, 0)
                        else np.array([1.0, 0.0])
                    )

                    patch_coords = CSSHoneycomb.__bulge_boundary_edge(
                        p_coords, gap_index, outward_dir, bulge, n_arc
                    )

                plt.fill(
                    patch_coords[:, 0],
                    patch_coords[:, 1],
                    alpha=0.3,
                    facecolor=p.colour,
                    edgecolor="black",
                    zorder=1
                )
            elif len(p) == 2:

                if p.missing_positions:
                    missing_mapped = coord_map(np.asarray(p.missing_positions, dtype=float))
                    missing_centre = missing_mapped.mean(axis=0)
                    p_centre = p_coords.mean(axis=0)
                    outward_dir = missing_centre - p_centre
                    norm = np.linalg.norm(outward_dir)
                    outward_dir = outward_dir / norm if not np.isclose(norm, 0) else np.array([1.0, 0.0])
                else:
                    outward_dir = np.array([1.0, 0.0])

                patch_coords = CSSHoneycomb.__bulge_edge_outward(
                    p_coords[0], p_coords[1], outward_dir, bulge, n_arc
                )

                plt.fill(
                    patch_coords[:, 0],
                    patch_coords[:, 1],
                    alpha=0.3,
                    facecolor=p.colour,
                    edgecolor="black",
                    zorder=1
                )

            elif len(p) == 1:
                centre = p_coords[0]

                # Estimate a sensible radius from the nearest other vertex
                distances = np.linalg.norm(mapped_coords - centre, axis=1)
                nonzero_distances = distances[distances > 1e-10]

                if len(nonzero_distances) > 0:
                    radius = 0.15 * np.min(nonzero_distances)
                else:
                    radius = 0.5

                circle = plt.Circle(
                    centre,
                    radius,
                    alpha=0.3,
                    facecolor=p.colour,
                    edgecolor="black",
                    zorder=1
                )

                plt.gca().add_patch(circle)


        for e in self.edges:

            e_coords = coord_map(
                np.asarray(
                    e.get_coordinates(),
                    dtype=float
                )
            )

            plt.plot(
                e_coords[:, 0],
                e_coords[:, 1],
                "-",
                c=CSSHoneycomb.PALETTE[e.colour],
                zorder=2
            )

        plt.scatter(
            mapped_coords[:, 0],
            mapped_coords[:, 1],
            c="black",
            zorder=3
        )

        plt.axis("equal")
        plt.gca().yaxis.set_inverted(True)
        
        if savepath is not None:
            plt.savefig(savepath)

        if show:
            plt.show()

    @staticmethod
    def __bulge_boundary_edge(p_coords, gap_index, outward_dir, bulge, n_arc=16):
        """
        Replace the edge at `gap_index` with a quadratic Bezier arc bulging
        towards `outward_dir`. `outward_dir` is derived per-plaquette from
        where its truncated (missing) vertex would have been, in whatever
        coordinate system is currently being plotted -- not from the
        polygon's own winding order, which is unreliable for boundary
        plaquettes that collapse to a near-zero-area sliver in square
        coordinates.
        """
        p_coords = np.asarray(p_coords, dtype=float)
        n = len(p_coords)

        if n < 3:
            raise ValueError("A polygon must contain at least 3 vertices.")

        i0 = gap_index
        i1 = (i0 + 1) % n
        p0, p1 = p_coords[i0], p_coords[i1]

        d = p1 - p0
        length = np.linalg.norm(d)
        if np.isclose(length, 0):
            return p_coords.copy()
        d_hat = d / length

        normal_1 = np.array([-d_hat[1], d_hat[0]])
        normal_2 = -normal_1

        outward = normal_2 if np.dot(normal_2, outward_dir) > np.dot(normal_1, outward_dir) else normal_1

        mid = (p0 + p1) / 2
        control = mid + outward * bulge * length

        t = np.linspace(0, 1, n_arc)[:, None]
        arc = (1 - t) ** 2 * p0 + 2 * (1 - t) * t * control + t ** 2 * p1

        remaining = [p_coords[(i1 + k) % n] for k in range(1, n - 1)]
        return np.vstack([arc, np.asarray(remaining)])
    
    @staticmethod
    def __bulge_edge_outward(p0, p1, outward_dir, bulge, n_arc=16):
        """
        For a weight-2 plaquette (two vertices, one edge), replace the
        straight edge with a single Bezier arc bulging outward. plt.fill
        closes the shape by drawing a straight line back from p1 to p0,
        so no return arc is needed.
        """
        p0 = np.asarray(p0, dtype=float)
        p1 = np.asarray(p1, dtype=float)

        d = p1 - p0
        length = np.linalg.norm(d)
        if np.isclose(length, 0):
            return np.vstack([p0, p1])
        d_hat = d / length

        normal_1 = np.array([-d_hat[1], d_hat[0]])
        normal_2 = -normal_1

        outward = normal_2 if np.dot(normal_2, outward_dir) > np.dot(normal_1, outward_dir) else normal_1

        control = (p0 + p1) / 2 + outward * bulge * length

        t = np.linspace(0, 1, n_arc)[:, None]
        arc = (1 - t) ** 2 * p0 + 2 * (1 - t) * t * control + t ** 2 * p1

        return arc
    