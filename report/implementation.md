We wrote our program in the Julia programming language, due to its design as a high-performance programming language designed for scientific computing. We also made a simple graphical interface (GUI) in Python that we use to generate interesting frameworks. We divided our code into the following 8 modules:
    * Beam1D.jl: Underlying solver for the single 1-dimensional case. Given a set of boundary conditions, this creates the system of matrices to solve, and has methods to solve this system in various ways.
    * single\_static.jl: Sets up the problem for the static single-beam case, and calls Beam1D.jl to solve it.
    * single\_dynamic\_newmark.jl: Sets up the problem for the dynamic single-beam case, and calls Beam1D.jl to solve it using the Newmark method.
    * single\_dynamic\_eigen.jl: Sets up the problem for the dynamic single-beam case, and calls Beam1D.jl to solve it using the eigenvalue method.
    * construct\textunderscore framework.py: GUI to create \emph{.ne} files, which we load into our framework solvers.
    * Beam2D.jl: Underlying solver for a framework of beams. Given a framework \emph{.ne} file, it sets up the system of matrices to solve.
    * framework\_static.jl: Sets up the problem for the static framework, and calls Beam2D.jl to solve it.
    * framework\_dynamic.jl: Sets up the problem for the dynamic framework, and calls Beam2D.jl to solve it. Can use either the Newmark or the Eigenvalue method.

Our ___Beam1D.jl___ file takes a set of problem parameters (boundary conditions, the grid for our solution, and functions describing the beam's stiffness and the force placed on the beam), and it calculates the system of the mass matrix $M_e$, the stiffness matrix $S_e$, and the constraint matrix $q_e$. It generates these by the finite element method, using cubic basis functions specified above and a 3 point Gaussian quadrature.
