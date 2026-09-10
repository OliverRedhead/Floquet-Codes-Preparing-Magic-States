# Floquet-Codes-Preparing-Magic-States

I can see our code working something like this:
1. A set of homology classes
    Given that this is a topological code it makes sense to make use of the mathematial tools we already have. I can see boundary and coboundary maps being helpful. I'm not sure what properties we will use other than their function of giving us the vertices (qubits) of an edge we want.
2. A set of qubit classes
    This was a really useful class in my heavy-hex, it helps us interact with stim which is often the hardest part of working with the library. I usually like to specify the utility of given qubits. For example qubits we initialise to $\ket0$ we might call Zdata, or ancillas on red edges we might call RedAncilla.
3. A set of helper methods (subroutines)
    This is usually the real important stuff. Things like parsing our grid to make sure we don't have any double qubits or methods to pick out all red edges to measure the relevant operators. The idea is that our homology and qubit classes are there to make these subroutines as easy as possible
4. A CSSHoneycomb class
    This is the class that will directly interact with stim. It helps me to have wrapper methods that do lots for us behind the scene so that we can interact with stim how we want to. This will be the real code specific methods, methods like measure_red_edges() which returns a string to be passed to stim (or possibly just a stim circuit).

    This is the really important class for debugging our actual circuit, making sure our measurement schedules are working properly. It will allow us to break our circuit into bits and should make it easier (a little bit easier, still going to be a pain)

# Progress
## 11/09/26
I have commited the start of this new project. Using basically the same qubit class as my previous projects and starting with a proper homology class. In my heavy-hex I used a similar concept in the Hexagon (hexagon.py) class. This basically was a homology class I just didnt know how what homology was. 

I can see the bulk of our code being easier with these tools in homology making it possible for us to think in terms of the schedule in the margaritta paper and let the code do all the hard work for us.

I can see some tweaks required for this class, but by in large I think it's structure makes sense. The pipeline for first part of stim initialisation looks something like this:
   1. Initialise the vertices we need, just use integers for `id`
   2. Initialise edges we need, this is easier said than done but we should only have to do this once and then forget about it
   3. Initialise plaquettes, again easier said than done but once these are set up we should be sweet
   4. We will need to come up with some coordinate system for our qubits so we can pair an edge with a qubit

**Possible adjustments/notes**
   - add a qubit attribute to our vertices to hold the relevant qubit
   - boundary conditions are a non-trivial problem at the moment, not sure if this method works well with them
   - I find it easier to start with the triangular superlattice and place cells given that, might use this to initialise everything. 
  
Still lots of work to do and lots of problems I havn't worked out how to solve, but this does give a solid starting point that should apply pretty generally.