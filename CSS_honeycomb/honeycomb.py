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
        assert isinstance(id, int | None)
    
    
        # qubit key must be integer in stim
        if not isinstance(id, int): 
            raise ValueError(f"id must be `int` not {type(id).__name__}")
        
        # qubit must have same key as this Vertex
        if isinstance(qubit, Qubit): 
            assert qubit.key == key, f"""`Vertex.key` and `Vertex.qubit.key` must match. 
                                        Cannot have vertex.key={id}, qubit.key={qubit.key}"""
            self.qubit = qubit

        elif qubit is None:
            self.qubit = Qubit(key, pos)
        
        else:
            raise ValueError(f"Vertex.qubit must be `Qubit` or `None`, not {type(qubit).__name__}")
        
        self.pos = pos

        super().__init__(id)

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

        super().__init__(id, v0, v1)
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
        vertices: Sequence[Vertex] | None = None
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

        super().__init__(id, edges)
        self.vertices = OneChain(vertices)

        if colour not in ["red", "green", "blue"]:
            raise ValueError(
                f"`colour` must be 'red', 'green', 'blue', not {colour}"
            )

        self.colour = colour

    def __str__(self):
        return f"Plaquette({self.key}) {self.colour}"

    def get_coordinates(
        self,
        ordered = False
    ) -> tuple[tuple[float, float], ...]:

        """Return coordinates of the unique vertices making up these edges.

        Parameters
        ----------
        order : bool = False
            If specified, order the vertices around their centre.
        """

        edges = self.boundary()

        vertices = {
            v
            for edge in edges
            for v in edge.boundary()
        }

        coords = np.array(
            [v.pos for v in vertices],
            dtype=float
        )

        if ordered is False:
            return tuple(map(tuple, coords))

        # Centre of the vertices
        centre = coords.mean(axis=0)

        # Angle of each vertex relative to the centre
        angles = np.arctan2(
            coords[:, 1] - centre[1],
            coords[:, 0] - centre[0]
        )

        # Increasing angle = anticlockwise
        indices = np.argsort(angles)
        coords = coords[indices]

        return tuple(map(tuple, coords))
    


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
        Colour choice is relatively simple in our square coordinate system.
        """

        plaquettes = []

        for i, v0 in enumerate(self.vertices):
            x, y = v0.pos

            if (x+y)%2 == 0 and x < self.ncols - 1 and y < self.nrows - 2:

                v1 = self.vertices[i + 1]
                v2 = self.vertices[i + self.ncols]
                v3 = self.vertices[i + self.ncols + 1]
                v4 = self.vertices[i + 2*self.ncols]
                v5 = self.vertices[i + 2*self.ncols + 1]

                plaquette = Plaquette(
                    len(plaquettes),
                    colour=Surface.__colour_plaquette(y),
                    vertices=[v0, v1, v2, v3, v4, v5]
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

    def plot_surface(self, coordinates="square"):

        """Visualise the surface."""

        if coordinates == "square":

            def coord_map(x):
                return np.asarray(x, dtype=float)

        elif coordinates == "hex":
            coord_map = self.square_to_hex # type: ignore

        else:
            raise ValueError(
                f"`coordinates` must be 'square' or 'hex', not {coordinates!r}"
            )

        # Plaquettes
        for p in self.plaquettes:

            p_coords = coord_map(
                np.asarray(p.get_coordinates(ordered=True), dtype=float)
            )

            plt.fill(
                p_coords[:, 0],
                p_coords[:, 1],
                alpha=0.3,
                facecolor=p.colour,
                edgecolor="black",
                zorder=1
            )

        # Edges
        for e in self.edges:

            e_coords = coord_map(
                np.asarray(e.get_coordinates(), dtype=float)
            )

            plt.plot(
                e_coords[:, 0],
                e_coords[:, 1],
                "-",
                c=Surface.PALETTE[e.colour],
                zorder=2
            )

        # Vertices
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