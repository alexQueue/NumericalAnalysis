# Description
## 1-dimensional case
The mathematical model of the static bending equation is given by
$$(EIw'')''(x) = q(x),\hspace{0.5cm}x\in\Omega$$
where the domain $\Omega = [0,L]$ is one-dimensional and $L$ is the length of the beam. Specifically we have
* $w(x)$: height of the neutral axis (bending curve) at $x$
* $E(x)$: Young's modulus
* $I(x)$: area moment of inertia
* $q(x)$: load at $x$
In addition,
* $M(x) = EIw''(x)$: bending moment at x
* $Q(x) = -(EIw'')'(x)$: shear force at x
can be used for setting the boundary conditions, influencing the bending curve.
## 2-dimensional case
In the 2-dimensional case each edge is modelled very similar to the 1-dimensional case, where we in local coordinates to the edge have a displacement in both the y-direction and x-direction where the displacement in the x-direction is modelled using the wave-equation
$\mu \ddot{v}- (EAv')' = f_1.$$
To get the displacement in global coordinates we simply rotate displacement according to the angle between the edge and the x-axis.

Since we have a network of _connected_ beams, we need some way to determine geometric constraints. Depending on what type each node is, we have different constraints. For all types we have that the end points between connecting beams should always be touching, and the angle between beams should be constant (i.e. stiffness of angle). For fixed nodes one of the end points needs to be fixed in space for all times.

These constraints are mathematically modelled as
* Stiffness of angle: $w_i'(\cdot,t) - w_j'(\cdot,t)=0$, where $\cdot$ is either the endpoint or startpoint of the edge, depending on the orientation between the two edges.
* Fixed bearing: $v_i(\cdot,t)=0$, $w_i(\cdot,t)=0$, $w_i'(\cdot,t)$.
* Movable bearing in direction $\vec{m}$: 
$$\left(\begin{bmatrix}\cos(\pi/2)&-\sin(\pi/2)\\\sin(\pi/2)&\cos(\pi/2)\end{bmatrix}\vec{m}\right)^T\begin{bmatrix}\cos(\phi_i)v_i(\cdot,t)-\sin(\phi_i)w_i(\cdot,t)\\ \sin(\phi_i)v_i(\cdot,t)+\cos(\phi_i)w_i(\cdot,t)\end{bmatrix}=0$$
* Linking condition:
$$\cos(\phi_i)v_i(\cdot,t)-\sin(\phi_i)w_i(\cdot,t) - \cos(\phi_j)v_j(\cdot,t)+\sin(\phi_j)w_j(\cdot,t)=0$$
$$\sin(\phi_i)v_i(\cdot,t)+\cos(\phi_i)w_i(\cdot,t) - \sin(\phi_j)v_j(\cdot,t)-\cos(\phi_j)w_j(\cdot,t)=0$$
