include("Beam1D.jl")
import CubicHermiteSpline
import LinearAlgebra
import Plots

#Initial condition
pars = (mu=x->1,EI=x->1,q=x->-1)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'M')=>0,
            (1,'Q')=>0)
grid = collect(LinRange(0,1,10))
grid_2 = collect(LinRange(0,1,40))

#Change load & solve dynamically
pars  = (mu=x->1, EI=x->1, q=x->-(x<0.5))
sys  = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
alpha = 2*[1, 1, 1, 1, 1, 1]#[4 , 20, 10, 50, 10]
beta = alpha

#compute eigenfunction and frequency numerically and analytically
eigenvals, eigenvec, w = Beam1D.get_numerical_eigenvalues(sys)
eigenfreq, eigenfunc_ana = Beam1D.get_analytic_eigenvalues(sys, pars, 1)

#Plot eigenfrequencies (eigenvalues)
Plots.plot(eigenvals[1:4], seriestype=:scatter, label = "numerical" )
Plots.plot!(eigenfreq[1:4], seriestype =:scatter, label = "analytical")
Plots.savefig("img/Eigenfrequencies_analytical_numerical.png")

#Plot analytical eigenfunctions
p1 = Plots.plot((eigenfunc_ana(grid_2))[1], xlim = [0, 1], title = "w1")
p2 = Plots.plot((eigenfunc_ana(grid_2))[2], xlim = [0, 1], title = "w2")
p3 = Plots.plot((eigenfunc_ana(grid_2))[3], xlim = [0, 1], title = "w3")
p4 = Plots.plot((eigenfunc_ana(grid_2))[4], xlim = [0, 1], title = "w4")
plt = Plots.plot(p1, p2, p3, p4, layout= (2, 2), legend= false)
Plots.savefig("img/Eigenfunctions_analytical.png")

#Plot numerical eigenfunctions
pp1 = Plots.plot(Beam1D.u_to_Vh(sys, eigenvec[:, 1]), title = "w1")
pp2 = Plots.plot(Beam1D.u_to_Vh(sys, eigenvec[:, 2]), title = "w2")
pp3 = Plots.plot(Beam1D.u_to_Vh(sys, eigenvec[:, 3]), title = "w3")
pp4 = Plots.plot(Beam1D.u_to_Vh(sys, eigenvec[:, 4]), title = "w4")
plt = Plots.plot(p1, p2, p3, p4, layout= (2, 2), legend= false)
Plots.savefig("img/Eigenfunctions_numerical.png")

plt = Plots.plot(p1, pp1, p2, pp2, p3, pp3,  p4, pp4, layout= (4, 2), legend= false)
Plots.savefig("img/Eigenfunctions_numerical_and_analytical.png")


n = 6
times = collect(LinRange(0,0.1,500))
eigenfunc_t = [(x,t) -> CubicHermiteSpline.CubicHermiteSplineInterpolation(
		            sys.problem.grid,eigenvec[1:2:end-2, i],eigenvec[2:2:end-2, i])(x) .* real(w(t, alpha[i], beta[i])) for i in 1:n]
     

anim = Plots.@animate for (j,t) in enumerate(times)
     #superposition
     res = 0
     for i in 1:n
          res = res .+ eigenfunc_t[i](grid_2, t)
     end
     Plots.plot(grid_2, res ,ylim=[-1.5,1.5], xlim =[0,1], label = "w(x)")
end

Plots.gif(anim, "img/beam_animation_eigenvals_superposition.gif", fps=15)  


for i in 1:n
     anim = Plots.@animate for (j,t) in enumerate(times)
          Plots.plot(grid_2, eigenfunc_t[i](grid_2, t) ,ylim=[-1.5,1.5], xlim =[0,1])
     end
     Plots.gif(anim, string("img/beam_animation_eigenvals_", i, ".gif"), fps=15) 
end 