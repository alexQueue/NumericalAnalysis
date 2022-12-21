include("Beam1D.jl")
import Plots
import LinearAlgebra
import LaTeXStrings


#Set system and parameters
#Set boundary conditions the same as assumed for the analytical case
pars = (mu=x->1,EI=x->1,q=x->10)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'H')=>0,
            (1,'G')=>0)
L = 1
grid = collect(LinRange(0,1,200))
prob = Beam1D.Problem(pars,BCs,grid)
sys  = Beam1D.System(prob)

IC   = [sys.Se\sys.qe zeros(sys.shape[2])]

#Set time discretization
ts = collect(LinRange(0,1,100))

#numerical solution
get_sol_num = Beam1D.solve_dy_eigen(sys)
sol_num     = get_sol_num(IC)

#analytical solution
freq_ana, ev_ana = get_eigen_cantilever(sys, pars, L, 4)
get_sol_ana = Beam1D.solve_dy_eigen_analytical(sys, freq_ana, ev_ana)
sol_ana = get_sol_ana(IC)

#plot 
anim = Plots.@animate for (j, t) in enumerate(ts)
    Plots.plot([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
    Plots.plot!(sol_num(t)...,0,1,ylim=[-0.2,0.2],color="black",label=false,linewidth=2)
    Plots.plot!(grid,sol_ana(t),ylim=[-0.2,0.2],color="grey",label=false,linewidth=2)
#Plots.savefig("presentation/gifs/superpos"*@sprintf("/frame%i.png", j))
end

Plots.gif(anim, "img/single/beam_animation_eigen_num_ana.gif", fps=15)

#numerical 
evals, evecs, freqs, modes = Beam1D.get_vibrations(sys)

# Plot analytical against numerical eigenfrequencies/i.e. eigemvalues
Plots.plot(freq_ana[1:4], seriestype=:scatter, label = "analytical", mc="blue")
Plots.plot!(freqs[1:4], seriestype=:scatter, label = "numerical", mc="red")
Plots.xlabel!(LaTeXStrings.L"j")#
Plots.ylabel!(LaTeXStrings.L" Eigenfrequency  $\omega_j$")
Plots.plot!(legend=:topleft)
Plots.savefig("img/single/Eigenfrequencies_analytical_numerical.svg")

