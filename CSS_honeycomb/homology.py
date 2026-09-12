"""
It seems really natural to me to have a series of container classes to help with the initialisation of our circuit.
Given that the CSS Honeycomb code is a topological code, it just makes sense to use these topological concepts to help
with the setup.

I plan to use boundary and coboundary maps to move between homologies, and if we need to it would make using the 
superlattice of the triangular toric code easier.
"""

from abc import ABC, abstractmethod


class nCell(ABC):
    """
    n-cell base class. We will use this base to define the vertex, edge and plaquette classes.

    `id` should be something hashable and unique within a dimension (e.g. a coordinate tuple, or a plain integer label). 
    Equality and hashing are based on id + type, so two cells built separately with the same id are treated as the same 
    cell -- important since chains are implemented as sets/frozensets of cells.
    """

    n = None

    def __init__(self, id) -> None:
        """
        Initialise nCell instance.

        Parameters
        ----------
        id : hashable
            should be something hashable and unique within a dimension. I imagine we will mostly be using unique integer
            identifiers.
        """
        self.id = id
        # initialise empty coboundary list
        self.coboundary_cells: set[nCell] = set()

    def __eq__(self, other) -> bool:
        return isinstance(other, type(self)) and self.id == other.id

    def __hash__(self) -> int:
        return hash((type(self).__name__, self.id))

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.id!r})"
    
    """
    and our boundary and coboundary maps as abstract methods
    """

    @abstractmethod
    def boundary(self) -> nChain | None:
        pass

    @abstractmethod
    def coboundary(self) -> nChain | None:
        pass
    

class ZeroCell(nCell):
    """
    In our case, a 0-cell is a vertex.
    """

    n = 0

    def __init__(self, id) -> None:
        """Initialise 0-cell instance."""
        super().__init__(id)
            
    def boundary(self) -> None:
        raise NotImplementedError("no -1-cells: boundary of a 0-cell is undefined here")
    
    def coboundary(self) -> OneChain:
        """
        vertex coboundary map

        Returns
        -------
        edge : Edge
            The edge 
        """
        return OneChain(self.coboundary_cells)

class OneCell(nCell):
    """
    In our case, 1-cell is an edge
    """

    n = 1
    
    def __init__(self, id, v0: ZeroCell, v1: ZeroCell) -> None:
        """
        Initialise 1-Cell instance (edge). An edge is defined by the vertices it connects.
        While initialising we set up the coboundary maps of the vertices

        Parameters
        ----------
        id : hashable
            unique identifier of this 1-cell. Should be consistent within all 1-chains
        v0 : ZeroCell
            A ZeroCell (vertex) that defines one end of this edge
        v1 : ZeroCell
            The other ZeroCell (vertex) that defines the other end of this edge
        """
        super().__init__(id)
        self.vertices = (v0, v1)
        v0.coboundary_cells.add(self)
        v1.coboundary_cells.add(self)

    def boundary(self) -> ZeroChain:
        """
        We should have set the boundary of this chain on initialisation.
        Make a ZeroChain out of that and return.
        """
        v0, v1 = self.vertices
        return ZeroChain({v0}) + ZeroChain({v1})
    
    def coboundary(self) -> nChain:
        """
        NOTE This will need to be setup when we construct the relevant 2-cell.
        """
        return TwoChain(self.coboundary_cells)
      
class TwoCell(nCell):
    """
    In our case, 2-cell is a plaquette
    """

    n = 2

    def __init__(self, id, edges: list[OneCell]) -> None:
        """
        Initialise TwoCell (plaquette) defined by its edges.

        Parameters
        ----------
        id : hashable
            Unique identifier for this plaquette, should be consistent will other plaquettes
        edges : list[OneCell]
            Edges that define this plaquette
        """
        super().__init__(id)
        self.edges = tuple(edges)
        for e in edges:
            e.coboundary_cells.add(self)

    def boundary(self) -> nChain:
        """return the edges that define this plaquette"""
        result = OneChain({}) # empty is good 
        for e in self.edges:
            result = result + OneChain({e})
        return result
    
    def coboundary(self) -> None:
        raise NotImplementedError("no 3-cells: coboundary of a 2-cell is undefined here")



class nChain(ABC):
    """
    n-cell base class. To apply our nice homological formalisms, we will need some base functions and things of our chains.
    Will add to this class as I need it.

    `frozenset` is a really useful structure here. Got this idea of Claude
    """

    n = None

    def __init__(self, cells) -> None:
        """
        Initialise nChain instance.

        Parameters
        ----------
        cells : iterable
            An iterable holding the nCells that make up the chain
        """
        self.cells = frozenset(cells)

    def __add__(self, other):
        if not isinstance(other, nChain):
            return NotImplemented
        return type(self)(self.cells.symmetric_difference(other.cells))
    
    def __radd__(self, other):
        """lets sum([chain1, chain2, ...]) work, since sum() starts at 0"""
        if other == 0:
            return self
        return NotImplemented

    def __eq__(self, other) -> bool:
        return isinstance(other, type(self)) and self.cells == other.cells
    
    def __hash__(self) -> int:
        return hash(self.cells)

    def __bool__(self) -> bool:
        return bool(self.cells)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({set(self.cells)})"
    
    def __iter__(self):
        return iter(self.cells)

    @abstractmethod
    def boundary(self):
        raise NotImplementedError

    @abstractmethod
    def coboundary(self):
        raise NotImplementedError
    

class ZeroChain(nChain):

    n = 0

    def boundary(self):
        raise NotImplementedError("no -1-cells: boundary of a 0-chain is undefined here")

    def coboundary(self) -> OneChain:
        result = OneChain({})
        for c in self.cells:
            result = result + c.coboundary()
        return result


class OneChain(nChain):

    n = 1

    def boundary(self) -> ZeroChain:
        result = ZeroChain({})
        for c in self.cells:
            result = result + c.boundary()
        return result

    def coboundary(self) -> TwoChain:
        result = TwoChain({})
        for c in self.cells:
            result = result + c.coboundary()
        return result

class TwoChain(nChain):

    n = 2

    def boundary(self) -> OneChain:
        result = OneChain({})
        for c in self.cells:
            result = result + c.boundary()
        return result

    def coboundary(self):
        raise NotImplementedError("no 3-cells: coboundary of a 2-chain is undefined here")



# --- #
# Check some basic properties
# --- #

if __name__ == "__main__":
    # six vertices around a hexagon
    verts = [ZeroCell(i) for i in range(6)]

    # six edges connecting them in a cycle
    edges = [OneCell(i, verts[i], verts[(i + 1) % 6]) for i in range(6)]

    # one face bounded by all six edges
    face = TwoCell("hex0", edges)

    face_chain = TwoChain({face})
    edge_boundary = face_chain.boundary()
    vertex_boundary = edge_boundary.boundary()

    print("boundary of face  :", edge_boundary)
    print("boundary of that  :", vertex_boundary)
    assert not vertex_boundary, "d^2 != 0 -- something's wired wrong!"
    print("d^2 = 0 check passed")

    # coboundary should be the dual statement: coboundary of a vertex's
    # coboundary should land back on faces containing it, etc.
    v0_coboundary = verts[0].coboundary()      # edges touching vertex 0
    print("coboundary of vertex 0:", v0_coboundary)
    print("coboundary of that    :", v0_coboundary.coboundary())