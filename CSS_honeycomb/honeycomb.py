"""
These classes will initalise our surface as a set of vertex, edge and plaquette objects 
and help us do some useful strings of operations
"""

import numpy as np
import matplotlib.pyplot as plt
from homology import ZeroCell, OneCell, TwoCell, ZeroChain, OneChain, TwoChain
from qubit import Qubit
from collections.abc import Sequence

"""
nCell subclasses to hold coordinates and colour information as well as inherit all the 
homology properties we want to use.
"""

class Vertex(ZeroCell):

    """
    - A vertex is a ZeroCell that also supports a data qubit.
    - A vertex has a unique key and a position. We will use integers as keys.
        It is important that keys are integers as a Vertex and its Qubit should have the same unique key.
    - A vertex has no colour.
    """

    def __init__(self, key: int, pos: tuple[float, float] = (0,0), *, qubit: Qubit | None = None) -> None:
        """
        Initialise a vertex instance.

        parameters
        ----------
        key : int 
            unique key of vertex, this will be used when defining stim circuit.
        pos : tuple[float, float] = (0, 0)
            position of vertex in our 2D space, top left is (0,0)
        qubit : Qubit | None = None
            choice to pass pre-initialised qubit. If qubit is `Qubit` then 
            the position parameter is ignored and vertex.qubit is set to 
            this qubit

        Note
        ----
        - If qubit is not specified, we just make one at the position we are interestesd
        - vertex key and qubit key **must** match, else we will have problems with stim.

        """
        assert isinstance(key, (int, type(None))), f"Vertex.key must be `int` or `None`. Not {type(key)}"
    
    
        # qubit key must be integer in stim
        if not isinstance(key, int): 
            raise ValueError(f"key must be `int` not {type(key).__name__}")
        
        # qubit must have same key as this Vertex
        if isinstance(qubit, Qubit): 
            assert qubit.key == key, f"""Vertex.key and Vertex.qubit.key must match. 
                                        Cannot have vertex.key={key}, qubit.key={qubit.key}"""
            self.qubit = qubit

        elif qubit is None:
            self.qubit = Qubit(key, pos)
        
        else:
            raise ValueError(f"Vertex.qubit must be `Qubit` or `None`, not {type(qubit).__name__}")
        
        self.pos = pos

        super().__init__(key)

    def get_coordinates(self):
        """return coordinates of this vertex as tuple"""
        return self.pos

class Edge(OneCell):

    """
    - An edge is a OneCell defined by the two Vertex instances that it joins.
    - We will use integers as unique keys, although we shouldn't use them (I think).
    - An edge does have a colour, which should be defined by the colour of plaquettes that it bridges.

    """

    def __init__(self, key, v0: Vertex, v1: Vertex, colour: str) -> None:
        if not isinstance(v0, Vertex) or not isinstance(v1, Vertex):
            raise ValueError(f"Edge must be handed two `Vertex` instances, not v0: {type(v0)} and v1: {type(v1)}")

        super().__init__(key, v0, v1)
        if colour not in ['red', 'green', 'blue']:
            raise ValueError(f"`colour` must be 'red', 'green', 'blue', not {colour}")
        self.colour = colour

    def get_coordinates(self):
        """return coordinates of the vertices that make up this edge as tuple"""
        v0, v1 = self.boundary()
        return (v0.pos, v1.pos) # type: ignore
    
    def get_indices(self):
        v0, v1 = self.boundary()
        return v0.key, v1.key

class Plaquette(TwoCell):

    """
    - A plaquette is a TwoCell defined by the Edge or Vertex instances that make up its boundary.
    - We will use integers as unique keys, although we shouldn't use them (I think)
    - A plaquetter does have a colour, we will use the convention that the top left hexagon is red,
        and the hexagon joining it to the lower right is green. This should uniquelly define the 
        colours of our surface.
    
    Note:
    This class allows us to initialise a plaquette given a list of its edges **or** its vertices (not both)
    In either case we compute the other too so a Plaquette instance holds a list of edges and vertices, not
    sure how useful this will be.
    """

    def __init__(
        self,
        key,
        colour: str,
        *,
        edges: Sequence[Edge] | None = None,
        vertices: Sequence[Vertex] | None = None,
        missing_positions: Sequence[tuple[float, float]] | None = None
    ) -> None:

        if vertices is not None: # if vertices are specified

            if edges is not None: # if edges are specified
                raise ValueError(
                    "You have specified both edges and vertices, "
                    "must only specify one or the other."
                )

            vertex_set = set(vertices)

            # Collect every edge incident to at least one supplied vertex
            candidate_edges = {
                edge
                for vertex in vertices
                for edge in vertex.coboundary()
                if isinstance(edge, Edge)
            }

            # Keep only edges whose two endpoints are both in the plaquette
            edges = tuple(
                edge
                for edge in candidate_edges
                if all(v in vertex_set for v in edge.vertices)
            )
        
        elif edges is not None: # if edges are specified

            vertex_set = {
                vertex
                for edge in edges
                for vertex in edge.boundary()
                if isinstance(vertex, Vertex)
            }

        else:
            raise ValueError(
                "cannot have edges and vertices None; "
                "must specify one or the other."
            )

        super().__init__(key, edges)
        self.vertices = OneChain(vertices)

        if colour not in ["red", "green", "blue"]:
            raise ValueError(
                f"`colour` must be 'red', 'green', 'blue', not {colour}"
            )

        self.colour = colour

        self.missing_positions = tuple(missing_positions) if missing_positions else ()

    def __str__(self):
        return f"Plaquette({self.key}) {self.colour}"
    
    def get_ordered_vertices(self) -> tuple["Vertex", ...]:
        """
        Return this plaquette's vertices, ordered anticlockwise around their
        centroid. Angles are computed from `.pos` (square coordinates), so
        this ordering is independent of whatever coordinate system is later
        used for plotting.
        """
        vertices = list(self.vertices)
        if not vertices:
            return tuple()

        coords = np.array([v.pos for v in vertices], dtype=float)
        centre = coords.mean(axis=0)
        angles = np.arctan2(coords[:, 1] - centre[1], coords[:, 0] - centre[0])
        order = np.argsort(angles)
        return tuple(vertices[i] for i in order)

    def get_coordinates(self, ordered=False) -> tuple[tuple[float, float], ...]:
        """Return coordinates of the unique vertices making up these edges."""
        if not ordered:
            vertices = {v for v in self.vertices}
            coords = np.array([v.pos for v in vertices], dtype=float)
            return tuple(map(tuple, coords))

        return tuple(v.pos for v in self.get_ordered_vertices())

    def get_boundary_gap_index(self) -> int | None:
        """
        Return the index i such that the edge between get_ordered_vertices()[i]
        and [i+1] (mod n) is *not* backed by a real Edge object -- i.e. it's
        the topological gap left by truncation at the surface boundary, which
        is what should be bulged outward. Returns None for a full, closed
        hexagon (no gap to speak of).
        """
        ordered = self.get_ordered_vertices()
        n = len(ordered)

        for i in range(n):
            v0, v1 = ordered[i], ordered[(i + 1) % n]
            if not (v0.coboundary_cells & v1.coboundary_cells):
                return i

        return None
    
    def __len__(self):
        """
        length of a plaquette is the number of vertices it supports
        """
        return len(self.vertices)


# ------------------------- #
#  Surface container class  #
# ------------------------- #

class Surface:
    """
    Owns the vertices, edges and plaquettes of a lattice, keyed by key.
    Responsible for construction, lookup, incidence queries and validation.
    """

    RGB = ['red', 'green', 'blue']

    """
    initialisation methods
    """

    def __init__(self, nrows : int = 4, ncols : int = 5) -> None:
        """
        Initialise Surface instance. Automatically generates vertices (and qubits) based on input nrows and ncols.
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
            colour = Surface.__colour_edge(v0.pos[1], v1.pos[1])
            return Edge(len(edges), v0=v0, v1=v1, colour=colour)

        edges = []
        for i, v0 in enumerate(self.vertices):
            x, y = v0.pos

            if y != 0:  # vertical edge, connects to the vertex above
                edges.append(make_edge(edges, self.vertices[i - self.ncols], v0))

            if (x + y) % 2 == 0 and x != self.ncols - 1:  # horizontal edge, checkerboard pattern
                edges.append(make_edge(edges, v0, self.vertices[i + 1]))
        return edges

    def __initialise_plaquettes(self) -> list[Plaquette]:
        """
        Initialise plaquettes, needs to be called after vertices are initialised.
        A plaquette is built starting from the qubit in the top left (closest to (0,0)).
        Anchors are allowed to sit one step outside the grid (x=-1 or y=-1) so that
        boundary plaquettes on the left/top are generated too, truncated to whichever
        of the 6 vertices actually exist. Weight-1 and weight-2 boundary stabilisers
        are kept.
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

                plaquette = Plaquette(
                    len(plaquettes),
                    colour=Surface.__colour_plaquette(y),
                    vertices=vertices,
                    missing_positions=missing_positions,
                )
                plaquettes.append(plaquette)

        return plaquettes

    """
    stim interaction methods
    """

    def initialise_qubits(self, coordinates='square') -> str:
        
        if coordinates != 'square':
            return self.__initialise_qubits_hex()
        
        string = ""
        for v in self.vertices:
            string += str(v.qubit) + "\n"
        
        return string
    
    def __initialise_qubits_hex(self):
        string = ""        
        hex_coords = Surface.square_to_hex(self.coords)
        for i, (x, y) in enumerate(hex_coords):
            string += f"QUBIT_COORDS({x}, {y}) {i}\n"
        return string

    def measure_edges(self, colour, flavour="Z"):
        """
        We can pick out all the edges of a certain colour and measure them easily
        """
        string = f"M{flavour} "
        target_edges = [e for e in self.edges if e.colour == colour]
        target_qubits = [e.get_indices() for e in target_edges]
        for t0, t1 in target_qubits:
            string += f"{t0} {t1} "
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
        return Surface.RGB[index]
    
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

                    patch_coords = Surface.__bulge_boundary_edge(
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

                patch_coords = Surface.__bulge_edge_outward(
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
                c=Surface.PALETTE[e.colour],
                zorder=2
            )

        mapped_coords = coord_map(
            np.asarray(self.coords, dtype=float)
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