include("Beam1D.jl")
import Plots

pars = (mu=x->2-x,EI=x->2-x,q=x->-100)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'H')=>0,
            (1,'G')=>0)
grid = collect(LinRange(0,1,4))

sys  = Beam1D.System(Beam1D.Problem(pars,BCs,grid))

IC   = [sys.Se\sys.qe zeros(sys.shape[2],2)]

modes, weights = Beam1D.solve_tr_num_eig(sys,IC)

Plots.plot(x->modes(x)'weights(0.),xlim=grid[[1,end]])

anim = Plots.@animate for t in LinRange(0,10,300)
     sol_t(x) = modes(x)'weights(t)
     Plots.plot(sol_t,ylim=[-1.5,1.5])
     Plots.plot!(sol_t,grid,seriestype=:scatter)
end

Plots.gif(anim, "beam_animation.gif", fps=15)