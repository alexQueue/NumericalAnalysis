include("Beam1D.jl")
import Plots

#Initial conditions
pars = (mu=x->1,EI=x->1,q=x->-10)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'M')=>0,
            (1,'Q')=>0)
grid = collect(LinRange(0,1,10))

prob = Beam1D.Problem(pars,BCs,grid)
sys  = Beam1D.System(prob)

IC   = [sys.Se\sys.qe zeros(sys.shape[2],2)]

#Change load & solve dynamically
prob.parameters.q(x) = 0
sys = Beam1D.System(prob)

times = collect(LinRange(0,10,50))
sols  = Beam1D.solve_dy_Newmark(sys,IC,times)

anim = Plots.@animate for (xs,ys) in sols
    Plots.plot(xs,ys,0,1,ylim=[-1.5,1.5],color="black",label=false,linewidth=2 )
    Plots.plot!([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
end

Plots.gif(anim, "img/beam_animation.gif", fps=15)

#TODO: Implement other configurations
