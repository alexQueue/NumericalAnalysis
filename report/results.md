# Results
## Analytical Solution
In this section we compare the analytical results of several load cases to our numerical results. For all cases a constant Youngs module $E(x) = E$ and a constant  area moment or inertia $I(x) = I$ is assumed. The numerical result is computed with 20 grid points and the length of the beam is $l = 1$. 
Fixing the beam on the left side leads to the boundary conditions 
    $$w(0) = 0 ,  \quad w'(x) = 0.
    \label{eq:Analy_BC_fixed_beam}$$
The first case is a simple constant load case, depicted in figure ref\{fig:Analytical_Sol_Beam_const_1\}, where we choose  $q_0 = 3$. By integrating the static bending equation ref\{eq:StaticBendingEq\} four times and using the boundary conditions ref\{eq:Analy_BC_fixed_beam\} and $Q(l) = 0$ and $M(l) = 0$  to find the integration constants 
we get the bending curve 
    $$EI w(x) = \frac{q_0 l^4}{24} \left( 6 \xi^2 - 4 \xi ^3 + \xi^4 \right)$$
where $\xi = \frac{x}{l}$. Looking at the bending curve $w(x)$ in figure ref\{fig:Analytical_Sol_Beam_const_2\}, the analytical result coincides with the numerical results and gives an $L_2$ error of TODO. 

FIGURE

The second case is displayed in figure ref\{fig:Analytical_Sol_Beam_decreas_1\}, where the load is decreasing with $q(x) = 3-3x$ and the same boundary conditions as in the first constant load example. The analytic result gives 
    $$EI w(x) = \frac{q_0 l^4}{120} \left( 10 \xi^2 -10\xi^3 + 5 \xi ^4 - \xi ^5 \right)$$
and is displayed in figure ref\{fig:Analytical_Sol_Beam_decreas_2\} together with the numerical result. 

## Validation
## Parameter study
