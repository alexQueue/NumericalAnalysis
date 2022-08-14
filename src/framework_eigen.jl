include("Beam2D.jl")
import Plots

problem = Beam2D.Problem("data/crane.ne")
sys = Beam2D.System(problem)

get_sol = Beam2D.solve_dy_eigen(sys)
u = sys.Se\sys.qe

IC = [u zeros(size(u)...,2)]
sol = get_sol(IC)
ts = collect(LinRange(0,100,100,))

undeformed = Beam2D.u_to_Vh(sys.problem,zeros(size(u)...))

anim = Plots.@animate for t in ts
    Plots.plot(undeformed[1],undeformed[2],0,1,color="black",label=false,linewidth=2,linestyle=:dot)
    Plots.plot!(sol(t)...,0,1,color="black",label=false,linewidth=2)
end

Plots.gif(anim, "img/framework_eigen.gif", fps=15)
