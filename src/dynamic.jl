include("Beam1D.jl")
import Plots

#Initial condition
pars = (mu=x->1,EI=x->1-x/2,q=x->-10)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'M')=>0,
            (1,'Q')=>0)
grid = collect(LinRange(0,1,10))

sys  = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
IC   = [sys.Se\sys.qe zeros(sys.shape[2],2)]

#Change load & solve dynamically
pars  = (mu=x->1, EI=x->1-x/2, q=x->-(x<0.5))
sys   = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
times = collect(LinRange(0,10,500))

sol = Beam1D.solve_tr_Newmark(sys,IC,times)

anim = Plots.@animate for (j,t) in enumerate(times)
    Plots.plot(sol[j],ylim=[-1.5,1.5])
    Plots.plot!(sol[j],sys.problem.grid,seriestype=:scatter)
end

Plots.gif(anim, "beam_animation.gif", fps=15)
