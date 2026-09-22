import numpy as np
import matplotlib.pyplot as plt
from qubit_containers import Vertex, Edge, Plaquette, BoundaryPlaquette

class CSSHoneycomb:
    """
    Owns the vertices, edges and plaquettes of a lattice, keyed by key.
    Responsible for construction, lookup, incidence queries and validation.
    """

    RGB = ['red', 'green', 'blue']

    """
    initialisation methods
    """

    def __init__(self, ncols : int = 5, nrows : int = 4) -> None:
        """
        Initialise CSSHoneycomb instance. Automatically generates vertices (and qubits) based on input nrows and ncols.
        NOTE that not all boundary conditions are accounted for. 

        Parameters
        ----------
        nrows : int = 4
            number of rows of qubits
        ncols : int = 5
            number of columns of qubits

        Notes
        -----
        - Coordinates are initialised in integer (square) system. This is more natural for the heavy hex. 
            I will try and make sure there are adequate methods to let us set up the circuit in both square 
            and hexagonal coordinates, but square makes most sense when initialising qubits and things.
        """

        self.nrows = nrows
        self.ncols = ncols

        self.coords, self.vertices = self.__initialise_vertices()
        self.edges = self.__initialise_edges()
        self.plaquettes = self.__initialise_plaquettes()

    def __initialise_vertices(self) -> tuple[list[tuple[int, int]], list[Vertex]]:
        """Initialise vertices, note that data qubits are initialised in this step too"""
        coords = [
            (x, y)
            for y in range(self.nrows) 
            for x in range(self.ncols)
            ]
        
        vertices = [
            Vertex(i, pos) 
            for i, pos in enumerate(coords)
            ]
        return coords, vertices

    def __initialise_edges(self) -> list[Edge]:
        """
        Initialise edges, needs to be called after vertices are initialised.
        An edge is defined by the two vertices it joins.
        """
        
        def make_edge(edges, v0, v1):
            colour = CSSHoneycomb.__colour_edge(v0.pos[1], v1.pos[1])
            ancilla_key = len(self.vertices) + len(edges)
            return Edge(len(edges), v0=v0, v1=v1, colour=colour, ancilla_key=ancilla_key)

        edges = []
        for i, v0 in enumerate(self.vertices):
            x, y = v0.pos

            if y != 0:  # vertical edge, connects to the vertex above
                edges.append(make_edge(edges, self.vertices[i - self.ncols], v0))

            if (x+y) % 2 == 0 and x != self.ncols - 1:  # horizontal edge, checkerboard pattern
                edges.append(make_edge(edges, v0, self.vertices[i + 1]))
        return edges

    def __initialise_plaquettes(self) -> list[Plaquette]:
        """
        Initialise plaquettes, needs to be called after vertices are initialised.
        A plaquette is built starting from the qubit in the top left (closest to (0,0)).
        Anchors are allowed to sit one step outside the grid (x=-1 or y=-1) so that
        boundary plaquettes on the left/top are generated too, truncated to whichever
        of the 6 vertices actually exist.
        """

        def get_vertex(x, y):
            if 0 <= x < self.ncols and 0 <= y < self.nrows:
                return self.vertices[y * self.ncols + x]
            return None

        offsets = [(0, 0), (1, 0), (0, 1), (1, 1), (0, 2), (1, 2)]

        plaquettes = []

        for y in range(-2, self.nrows):
            for x in range(-1, self.ncols):
                if (x + y) % 2 != 0:
                    continue

                candidates = [get_vertex(x + dx, y + dy) for dx, dy in offsets]
                vertices = [v for v in candidates if v is not None]
                missing_positions = [
                    (x + dx, y + dy)
                    for (dx, dy), v in zip(offsets, candidates)
                    if v is None
                ]

                if not vertices:
                    continue

                colour = CSSHoneycomb.__colour_plaquette(y)

                if missing_positions:
                    plaquette = BoundaryPlaquette(
                        len(plaquettes),
                        colour=colour,
                        vertices=vertices,
                        missing_positions=missing_positions,
                    )
                else:
                    plaquette = Plaquette(
                        len(plaquettes),
                        colour=colour,
                        vertices=vertices,
                    )

                plaquettes.append(plaquette)

        return plaquettes

    """
    stim interaction methods
    """

    # circuit initialisation

    def initialise_circuit(self, coordinates='square') -> str:
        """
        Converts our CSSHoneycomb into stim language.

        NOTE For now, the circuit only works in square coordinates.
        """
        
        if coordinates != 'square':
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
        hex_coords = CSSHoneycomb.square_to_hex(self.coords)
        for i, (x, y) in enumerate(hex_coords):
            string += f"QUBIT_COORDS({x}, {y}) {i}\n"
        return string

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

        if flavour == "Z":
            # ZZ: ancilla -> data
            return ancilla.key, vertex.key

        elif flavour == "X":
            # XX: data -> ancilla
            return vertex.key, ancilla.key

        else:
            raise ValueError("flavour must be 'X' or 'Z'")


    def measure_edges(self, colour, flavour="Z"):
        """
        We can pick out all the edges of a certain colour and measure them easily.
        """
        string = ""

        target_edges = [e for e in self.edges if e.colour == colour]
        ancilla_keys = [e.ancilla.key for e in target_edges]

        # Initialise ancillas in |0>
        string += f"# prepare ancillas in {flavour} basis\nR "
        for a in ancilla_keys:
            string += str(a) + " "
        string += "\nTICK\n"
        
        # If measuring X stabilizer, put ancillas X basis
        if flavour == "X":
            string += "H "      # Hadamard
            for a in ancilla_keys:
                string += str(a) + " "
            string += "\nTICK\n"

        # First round of CNOTs between v0 and ancilla
        string += "# fold stabilizers\nCNOT "
        for edge in target_edges:
            c, t = CSSHoneycomb.__get_CNOTs(edge, flavour=flavour, r=0)
            string += f"{c} {t} "
        string += "\nTICK\n"

        # Second round of CNOTs between v1 and ancilla
        string += "CNOT "
        for edge in target_edges:
            c, t = CSSHoneycomb.__get_CNOTs(edge, flavour=flavour, r=1)
            string += f"{c} {t} "
        string += "\nTICK\n"

        # If measuring X stabilizer, put ancillas back in Z basis
        if flavour == "X":
            string += "# rotate ancillas back to Z basis\nH "
            for a in ancilla_keys:
                string += str(a) + " "
            string += "\nTICK\n"

        # measure ancillas
        string += "M "
        for a in ancilla_keys:
            string += str(a) + " "
    
        return string

    # state preparation TODO

    def prepare_qubits(self) -> str:
        """
        This method should prepare our states states following our injection protocol

        NOTE For now, we initialise everying in |0>
        """
        qubit_keys = [v.key for v in self.vertices]
        
        string = "# initialise qubits in |0>\nR "
        for key in qubit_keys:
            string += str(key) + " "
        return string



    """
    helper methods
    """

    @staticmethod
    def __colour_plaquette(y):
        """
        for colouring plaquettes. Very simple condition in our choice of square coordinates
        """
        index = int(y%3)
        return CSSHoneycomb.RGB[index]
    
    @staticmethod
    def __colour_edge(y0, y1):
        """
        For colouring edges. Very simple condition in our choice of square coordinates.
        """
        rgb = ['red', 'green', 'blue']
        if y0 == y1:
            return rgb[int((y1 - 1) % 3)]
        return rgb[int( (y0 + 1) % 3 )]

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

        even = (y % 2 == 0)

        x[even] += 2 * a * (x[even] // 2)

        x[~even] += (
            a * ((x[~even] + 1) // 2)
            + a * ((x[~even] - 1) // 2)
        )

        y *= 2 * a

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

    def plot_surface(self, coordinates="square", bulge=0.35, n_arc=16):

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

            mapped_coords = coord_map(
               np.asarray(self.coords, dtype=float)
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