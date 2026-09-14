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

## 12/9/26
I have done a bit of work on our Surface class today, we can now set up a circuit! I am trying to implement methods to visualise things as I go and make sure everything is working properly. The main things I have done today
   - Vertex, edge and plaquette classes to hold colours and coordinates with some helpful getter methods that really will help with our operations
   - Initialisation methods in the Surface class to set up a hexagonal surface. I am using 'square coordinates' to set things up (they are much easier to find patterns in and allow us to use binary operations) and then I have written a map method which takes us from the square coords to the hexagonal coords. It may be useful to think in both systems so it will be nice to be able to move between them.
   - Once the homology was set up it was pretty easy to add some visualisation methods which draw our surface using matplotlib, this may be helpful for debugging and making sure things are where they are supposed to be. 
   - I started making sure we can use our homology properties to easily apply operations. We can! It is really easy to measure the qubits on red edges for example. Hopefully we can automate the application of gates and things too. 

**Some things to implement and worry about:**
   - Detectors and operation classes. Last year I had a set of classes for the operations. I'm not sure if we want to do the same thing, or just keep track of things a different way. I found that because of the way that detectors are indexed in stim they are some of the most difficult operations to work with (and some of the most important!). Last year I overcame this by just holding onto a record of all measurements in the circuit and having to call on that every time I wanted to place a detector.
   - Boundary Conditions. Honestly in writing this I havn't looked much at how we consider our boundaries, but it it obviously something we need to work out and something I havn't considered when writing this.
   - Qubit indices in visualisation methods. Should be able to label vertices by their id to help with debugging
   - Ancilla qubits. For heavy hex (not sure about phenomenological model) we will have ancillas on the "heavy" part of our hexagons. Then we apply cnots to flow stabilizers into measurements. This should actually be relatively easy with our homology tools!
# 14/9/26
   - fixed an error in the Plaquette initialisation that initialised it with vertices rather than edges. Thought it would be better to restructure some of the other classes like TwoCell. Can now initialise a Plaquette given its edges or vertices.
  
Also started work on operations classes which make our life a little easier when applying detectors and things in our circuit. In the process we have lost the ability to use our circuit in the hexagonal (non-integer) coordinates. If we really want our circuit to work in both systems, we should just be able to set up coordinate maps but will worry about that later. Just note that using our operations on hexagonal coordinates don't work at the moment.
