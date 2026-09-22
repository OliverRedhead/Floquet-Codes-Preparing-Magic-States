"""
These classes will initalise our surface as a set of vertex, edge and plaquette objects 
and help us do some useful strings of operations
"""

import numpy as np
from homology import ZeroCell, OneCell, TwoCell, ZeroChain, OneChain, TwoChain
from qubit import Qubit, Ancilla
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
    - Edges have an ancilla at the middle of the edge
    """

    def __init__(self, key, v0: Vertex, v1: Vertex, colour: str, ancilla_key: int) -> None:
        """
        Notes
        -----
        - In creating the ancilla on this edge, we just use the key passed for this edge, this may cause issues.
        """
        if not isinstance(v0, Vertex) or not isinstance(v1, Vertex):
            raise ValueError(f"Edge must be handed two `Vertex` instances, not v0: {type(v0)} and v1: {type(v1)}")

        super().__init__(key, v0, v1)
        if colour not in ['red', 'green', 'blue']:
            raise ValueError(f"`colour` must be 'red', 'green', 'blue', not {colour}")
        self.colour = colour
        
        # create ancilla
        x0, y0 = v0.pos
        x1, y1 = v1.pos
        xa = (x0 + x1) / 2
        ya = (y0 + y1) / 2
        self.ancilla = Ancilla(key=ancilla_key, pos=(xa, ya))
        

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
    - A plaquette does have a colour, we will use the convention that the top left hexagon is red,
        and the hexagon joining it to the lower right is green. This should uniquely define the
        colours of our surface.

    Note:
    This class allows us to initialise a plaquette given a list of its edges **or** its vertices (not both)
    In either case we compute the other too so a Plaquette instance holds a list of edges and vertices, not
    sure how useful this will be.

    A plaquette that is truncated by the surface boundary (i.e. is missing vertices or edges)
    should be represented by `BoundaryPlaquette` instead.
    """

    def __init__(
        self,
        key,
        colour: str,
        *,
        edges: Sequence[Edge] | None = None,
        vertices: Sequence[Vertex] | None = None,
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

            # `vertices` is None on this path, so fill it in from the edges
            vertices = tuple(vertex_set)

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
            vertices = {v for v in self.vertices if isinstance(v, Vertex)}
            coords = np.array([v.get_coordinates() for v in vertices], dtype=float)
            return tuple(map(tuple, coords))

        return tuple(v.pos for v in self.get_ordered_vertices())

    def __len__(self):
        """
        length of a plaquette is the number of vertices it supports
        """
        return len(self.vertices)

class BoundaryPlaquette(Plaquette):
    """
    A special type of plaquette that can contain less than 6 vertices and edges.

    These plaquettes lie on the boundary (hence BoundaryPlaquette).
    We need some special cases for these to make sure boundary conditions are satisfied.

    Attributes
    ----------
    missing_positions:
        Positions of the vertices that would complete the hexagon but were
        truncated away by the surface boundary.
    """

    def __init__(
        self,
        key,
        colour: str,
        *,
        edges: Sequence[Edge] | None = None,
        vertices: Sequence[Vertex] | None = None,
        missing_positions: Sequence[tuple[float, float]] | None = None,
        flavour: str | None = None
    ) -> None:
        
        super().__init__(key, colour, edges=edges, vertices=vertices)

        self.missing_positions = tuple(missing_positions) if missing_positions else ()
        self.flavour = flavour

    def __str__(self):
        return f"BoundaryPlaquette({self.key}) {self.colour}"

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