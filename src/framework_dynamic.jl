include("Beam2D.jl")
import Plots

problem = Beam2D.Problem("data/crane.ne")
sys  = Beam2D.System(problem)

times = collect(LinRange(0,100,100))

u = sys.Se\sys.qe
IC = [u zeros(size(u)...,2)]

sys.qe[:] .= 0

xy = Beam2D.solve_dy_Newmark(sys,IC,times)
undeformed = Beam2D.u_to_Vh(sys.problem,zeros(size(u)...))

anim = Plots.@animate for (j,t) in enumerate(times)
    Plots.plot(undeformed[1],undeformed[2],0,1,color="black",label=false,linewidth=2,linestyle=:dot)
    Plots.plot!(xy[j][1],xy[j][2],0,1,color="black",label=false,linewidth=2)
end
gif(anim, "img/framework_dynamic.gif", fps=fps)
