include("Beam1D.jl")
include("Diagnostics.jl")
import Plots

#Initial condition
pars = (mu=x->1,EI=x->1,q=x->10-(x<0.5))
BCs  = Dict((0,'H')=>0,
            (0,'G')=>1,
            (1,'M')=>2,
            (1,'Q')=>3)
n = 10
grid = collect(LinRange(0,1,n))
h = grid[2] - grid[1]

times = collect(LinRange(0,10,500))

sys  = Beam1D.build(Beam1D.Problem(pars,BCs,grid))
IC   = [sys.E\sys.f zeros(sys.shape[2],2)]

#Change load & solve dynamically
pars  = (mu=x->1, EI=x->1, q=x->-(x<0.5))
sys   = Beam1D.build(Beam1D.Problem(pars,BCs,grid))

sol, ut = Beam1D.solve_tr(sys,IC,times)

anim = Plots.@animate for (j,t) in enumerate(times)
    Plots.plot(sol[j],ylim=[-1.5,1.5])
    Plots.plot!(sol[j],sys.problem.grid,seriestype=:scatter)
end

Plots.gif(anim, "img/beam_animation.gif", fps=15)

Diagnostics.print_diagnostics(ut[:, 1], h)
