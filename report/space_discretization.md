# Spacial Discretization
To solve the static bending equation with the Finite Element Method we first find the weak formulation of the problem, then apply the Galerkin reduction, insert the boundary conditions and solve the resulting linear equation system.

Our trial and function space ``V`` includes piecewise twice differentiable functions ``V:=C^{2,p}(\Omega) = \{\phi:\Omega\to\mathbb{R}|\phi\in C(\Omega), \phi'\in C(\Omega), \phi''\in C^p(\Omega)\}``, since we need to consider the fourth order derivative of w(x). We begin by intergrating the static bending equation over the domain ``\Omega`` and multiplying with the test function ``\psi(x)\in V``

```math
\begin{align}
    \int _\Omega (EIw'')''(x) \psi (x) dx &= \int _\Omega q(x) \psi (x) dx, \quad \forall \psi \in V
\end{align}
```
Applying partial integration twice then yields
```math
\begin{align*}
   \int _\Omega (EIw'')'' \psi dx &= - \int_\Omega (EIw'')'\psi' dx + [(EIw'')' \psi]_\Omega \\
      & = \int _\Omega EIw'' \psi '' dx + \left[ (EIw'')'\psi \right]_\Omega - \left[ EIw''\psi' \right ]_\Omega \\
         & = \int _\Omega EIw'' \psi '' dx  \underbrace{Q(L) \psi(L) - Q(0)\psi(0) + M(L)\psi(L)-M(0)\psi(0) }_{b(\psi)} \\
            \Longleftrightarrow \quad \underbrace{ \int _\Omega EIw'' \psi '' dx }_{a(w, \psi)}  & =  \underbrace{\int_\Omega q \psi dx - b(\psi) }_{F(\psi)}
\end{align*}
```
```math
\begin{align*}
   \int _\Omega (EIw'')'' \psi dx &= - \int_\Omega (EIw'')'\psi' dx + [(EIw'')' \psi]_\Omega \\
   & = \int _\Omega EIw'' \psi '' dx + \left[ (EIw'')'\psi \right]_\Omega - \left[ EIw''\psi' \right]_\Omega \\
   & = \int _\Omega EIw'' \psi '' dx  \underbrace{Q(L) \psi(L) - Q(0)\psi(0) + M(L)\psi(L)-M(0)\psi(0) }_{b(\psi)} \\
   \Longleftrightarrow \quad \underbrace{ \int _\Omega EIw'' \psi '' dx }_{a(w, \psi)}  & =  \underbrace{\int_\Omega q \psi dx - b(\psi) }_{F(\psi)}
\end{align*}
```
Now we have the weak formulation to to find ``w``, s.t. ``a(w, \psi) = F(\psi) \quad \forall \psi \in V``.
To reduce the problem to a finite dimension we use the Galerkin substitution of ``V`` for ``V_n``, spanned by ``n`` basis functions ``\phi_i``, so that ``w`` becomes ``w_{n}(x) = \sum_{k=1}^{n} \alpha_k \phi_k(x)`` and ``\psi_j`` becomes ``\phi_j``, ``\forall \quad j \in \{1, ..., n \}``:
```math
\begin{equation*}
    \sum_k \alpha_k \left( \int_{\Omega} ( EI\phi_k'' \phi_j'') dx \right) = \int_{\Omega} q  \phi_j dx + b(\psi) \quad \forall j, \quad j, k \in \{1, ..., n\}
\end{equation*}
```
Not sure about this part since we have cubic basis functions: 
The basis functions are cubic on the elements ``\Omega_k = [x_k, x_{k+1})``., choosen s.t. ``\phi_i(x_j) = \delta_{ij}`` (not sure here). Hence on each element we have
```math
\begin{align*}
    \sum_k \alpha_k \int_{\Omega_{i}} ( EI\phi_k'' \phi_j'') dx  = \sum_k \int_{\Omega_{i}} q(x) \phi_j dx + b(\psi) \quad \forall j, \quad j, k \in \{1, ..., n\} \cap \{k, k+1, k?\}
\end{align*}
```
Each element of the stiffness matrix ``\mathbf{S}`` is calculated by
```math
\begin{align*}
    s_{jk} = \int_{0}^L E(x)I(x)\phi''_{k}\phi''_{j}.
\end{align*}
```

In the static case, the resulting system 
```math
\begin{equation*}
    \mathbf{S_e} \mathbf{u} =  \mathbf{q_e}
\end{equation*}
```
is solved and the solution passed to the cubic basis functions, which are used for the ansatz space ``V_h``. Using piecewise polynomials of degree 3 is necessary, since the static bending equation is of order 4, s.t. the weak formulation requires a twice differentiable basis function to compute the stiffness matrix. 
The ansatz space ``V_h`` with ``h = \frac{L}{n-1}``, ``x_i= h(i-1), i= 1, ..,`` and ``n\leq 2`` is defined as
```math
\begin{equation*}
    V_h =\{ \phi \in V | \phi_{]x_i, x_{i+1}[} , i= 1, .., n-1\}
\end{equation*}
```
where each function ``\phi_i \in V_h`` is uniquely defined by ``u_{2i-1} = \phi(x_i)`` and ``u_{2i} = \phi'(x_i)``. As a consequence each function can be determined by the weighted linear combination of the basis functions  ``\phi_i \in V_h`` 
```math
\begin{align*}
    \phi = \sum_{k=1}^{2n} u_k \phi_k = \sum_{k=1}^{n} (u_{2i-1}\phi_{2i-1} + u_{2i} \phi_{2i}).
\end{align*}
```
The basis functions ``\phi_i`` with ``i=2,.., n-1 `` are defined as follows
```math
\begin{align*}
    \phi_{2i-1}(x) = \begin{cases} 
    \Bar{\phi_3}(\frac{x-x_{i-1}}{h})  \quad & x \in [x_{i-1}, x_i] \\
    \Bar{\phi_1}(\frac{x-x_i}{h}) \quad & x \in [x_i, x_{i+1}] \\
    0 \quad  &\text{otherwise}
    \end{cases}
    \quad
    \phi_{2i}(x) = \begin{cases} 
    h \Bar{\phi_4}(\frac{x-x_{i-1}}{h})  \quad & x \in [x_{i-1}, x_i] \\
    h \Bar{\phi_2}(\frac{x-x_i}{h}) \quad & x \in [x_i, x_{i+1}] \\
    0 \quad  &\text{otherwise.}
    \end{cases}
\end{align*}
```
At the boundary elements, we get
```math
\begin{align*}
    \phi_{1}(x) = \begin{cases} 
    \Bar{\phi_1}(\frac{x}{h})  \quad & x \in [0, h] \\
    0 \quad  &\text{otherwise}
    \end{cases}
    \quad
    \phi_{2}(x) = \begin{cases} 
    h \Bar{\phi_2}(\frac{x}{h})  \quad & x \in [0, h] \\
    0 \quad  &\text{otherwise}
    \end{cases}
\end{align*}
```
and 
```math
\begin{align*}
    \phi_{2n-1}(x) = \begin{cases} 
    \Bar{\phi_3}(\frac{x-x_{n-1}}{h})  \quad & x \in [x_{n-1}, L] \\
    0 \quad  &\text{otherwise}
    \end{cases}
    \quad
    \phi_{2n}(x) = \begin{cases} 
    h \Bar{\phi_4}(\frac{x-x_{n-1}}{h})  \quad & x \in [x_{n-1}, L] \\
    0 \quad  &\text{otherwise.}
    \end{cases}
\end{align*}
```
The four form functions ``\Bar{\phi_i}, i = 1, ..,4`` used for the  basis functions are polynomial functions of degree 3 and determined by 
```math
\begin{align}
    \Bar{\phi_1}(\xi) &= 1-3\xi^2+ 2\xi^3 \\
    \Bar{\phi_2}(\xi) &= \xi(\xi-1)^2 \\
    \Bar{\phi_3}(\xi) &= 3\xi^2- 2\xi^3 \\
    \Bar{\phi_4}(\xi) &= \xi^2(\xi-1). \\
\end{align}
```
