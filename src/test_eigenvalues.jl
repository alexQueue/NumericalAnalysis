include("Beam1D.jl")
import Plots

pars = (mu=x->1,EI=x->1,q=x->-10)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'M')=>0,
            (1,'Q')=>0)
grid = collect(LinRange(0,1,5))

sys  = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
IC   = [sys.Se\sys.qe zeros(sys.shape[2])]

modes, weights = Beam1D.solve_tr_num_eig(sys)

weights_IC = weights(IC)

anim = Plots.@animate for t in LinRange(0,10,300)
     Plots.plot(x->modes(x)'weights_IC(t),ylim=[-1.5,1.5])
end

Plots.gif(anim, "beam_animation.gif", fps=15)