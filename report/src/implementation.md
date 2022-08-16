# Implementation
We wrote our program in the Julia programming language, due to its design as a high-performance programming language designed for scientific computing. We also made a simple graphical interface (GUI) in Python that we use to generate interesting frameworks. We divided our code into the following 9 modules:
* Beam1D.jl: Underlying solver for the single 1-dimensional case. Given a set of boundary conditions, this creates the system of matrices to solve, and has methods to solve this system in various ways.
* single\_static.jl: Sets up the problem for the static single-beam case, and calls Beam1D.jl to solve it.
* single\_dynamic.jl: Sets up the problem for the dynamic single-beam case, and calls Beam1D.jl to solve it using the Newmark method.
* single\_eigen.jl: Sets up the problem for the dynamic single-beam case, and calls Beam1D.jl to solve it using the eigenvalue method.
* construct\_framework.py: GUI to create \emph{.ne} files, which we load into our framework solvers.
* visualize\_framework.jl: Takes one of the generated .ne files to visualize it with each node colored corresponding to the type of the node.
* Beam2D.jl: Underlying solver for a framework of beams. Given a framework \emph{.ne} file, it sets up the system of matrices to solve.
* framework\_static.jl: Sets up the problem for the static framework, and calls Beam2D.jl to solve it.
* framework\_dynamic.jl: Sets up the problem for the dynamic framework, and calls Beam2D.jl to solve it using the Newmark method.
* framework\_eigen.jl: Sets up the problem for the dynamic framework, and calls Beam2D.jl's solve\_dy\_eigen to solve it. 
## Detailed descriptions
### BeamXD.jl (X=1 or 2)
The _BeamXD.jl_ files both make use of a _Problem_ struct and a _System_ struct. The _Problem_ struct holds all the information about the problem statement, whereas the _System_ struct defines the matrices that are used for solving the problem.

The _System_ struct is similar in both cases as it calculates the matrices $M_e$, $S_e$ and $q_e$ using the methods described in the math parts. (TODO: Reformulate that).

For plotting we need to convert the resulting displacement vector $u$ into parametric curves. We chose to use --- curves. The function _u\_to\_Vh_ does that and returns a vector/scalar of the parametric curves.

The functions _solve\_dy\_newmark_ and _solve\_dy\_eigen_ returns a function to calculate the solution of the problem at a specific time step using the newmark and eigenvalue method respectively.

### Beam1D.jl
Our _Beam1D.jl_ file takes a set of problem parameters (boundary conditions, the grid for our solution, and functions describing the beam's stiffness and the force placed on the beam), and it calculates the system of the mass matrix $M_e$, the stiffness matrix $S_e$, and the constraint matrix $q_e$. It generates these by the finite element method, using cubic basis functions specified above and a 3 point Gaussian quadrature.

### Beam2D.jl
The _Beam2D.jl_ file instead takes in a .ne (node-edge) file that specifies the problem structure using coordinates of nodes and node connections, as well as the type of each node (Free/Fixed/Force/Movable). The file format is 

```
NODES
<x-coord> <y-coord>
<x-coord> <y-coord>
...
EDGES
<node 1> <node 2> [optional: params=[<fnc>,<fnc>,<fnc>,<fnc>],gp=x]
<node 1> <node 2>
...
TYPES
1 <FIXED/FREE/FORCE/MOVABLE>
2 <FIXED/FREE/FORCE/MOVABLE>
3 <FIXED/FREE/FORCE/MOVABLE>
...
PARAMETERS
E <fnc>
I <fnc>
A <fnc>
mu <fnc>
```
where the optional ´´params´´ defines specific material properties for that edge, and the optional ´´gp´´ defines a specific amount of gridpoints at that edge.

Using a file in this format the function  _Problem_ creates a new instance of the struct. This struct has three members (TODO:fix this word), _nodes_, _edges_, and the shape of the problem. The _nodes_ and _edges_ members are vectors of the structs _Node_ and _Edge_, which holds information about the specific edge and node. 
