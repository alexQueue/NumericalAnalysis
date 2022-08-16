# Time discretization
To solve the dynamic case of the beam bending equation the Newmark Method was used. This method finds an approximation of the solution to the equation
$$ f(\ddot{u}(t),\dot{u}(t),u(t),t)=0.$$
In the case considered in this project $f(\ddot{u},\dot{u},u,t) = M\ddot{u}+D\dot{u}+Su-q(t)$.
Considering a partition of the interval we are trying to find the solution in, $t_1\le t_2\le \dots \le t_n$, such that $u(t_1) = u_1$, $\dot{u}(t_1) = \dot{u}_1$ and $\ddot{u}(t_1) = \ddot{u}_1$. 
Then the algorithm is as follows:
* For $j=1$ to $n-1$:
** Compute 
$$u_j^* = u_j + \dot{u}_jh_j + \left(\frac{1}{2}-\beta\right)\ddot{u}_jh_j^2$$
and
$$\dot{u}_j^* = \dot{u}_j + \left(1-\gamma\right)\ddot{u}_jh_j^2$$
where $h_j = t_{j+1}-t_j$.
<!-- $$u_{j+1}=u_j+\dot{u}_jh_j + \left(\left(\frac{1}{2}-\beta\right)\ddot{u}_j + \beta\ddot{u}_{j+1}\right)h_j^2$$ --> 
** Solve 
$$\left(M+\gamma h_jD+\beta h_j^2S)\ddot{u}_{j+1} = q(t_{j+1})-D\dot{u}_j^*-Su_j^*$$
** Compute the new iterations 
$$u_{j+1} = u_j^* +\beta \ddot{u}_{j+1}h_j^2$$
and
$$\dot{u}_{j+1} = \dot{u}_j^* + \gamma \ddot{u}_{j+1}h_j.$$
