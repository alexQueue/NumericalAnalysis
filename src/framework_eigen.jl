include("Beam2D.jl")
import Plots

problem = Beam2D.Problem("data/framework3.ne")
sys = Beam2D.System(problem)

get_sol = Beam2D.solve_dy_eigen(sys)
u = sys.Se\sys.qe

IC = [u zeros(size(u)...,2)]
sol = get_sol(IC)
ts = collect(LinRange(0,100,100,))

anim = Plots.@animate for t in ts
    Plots.plot(sol(t)...,0,1,color="black",label=false)
end

Plots.gif(anim, "img/framework_eigen.gif")
