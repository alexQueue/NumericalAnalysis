include("Beam2D.jl")
import Plots
using Printf

problem = Beam2D.Problem("data/crane.ne")
sys = Beam2D.System(problem)

get_sol = Beam2D.solve_dy_eigen(sys)
u = sys.Se\sys.qe

IC = [u zeros(size(u)...,2)]
sol_eigen = get_sol(IC)
ts = collect(LinRange(0,100,500,))

sys.qe[:] .= 0
xy = Beam2D.solve_dy_Newmark(sys,IC,ts)

undeformed = Beam2D.u_to_Vh(sys.problem,zeros(size(u)...))

newmark_label = hcat("Newmark",fill("",1,99))
eigen_label = hcat("Eigen",fill("",1,99))
anim = Plots.@animate for (j,t) in enumerate(ts)
    Plots.plot(undeformed[1],undeformed[2],0,1,color="black",label=false,legend=:topright,
        linewidth=2,linestyle=:dot)
    Plots.plot!(xy[j][1],xy[j][2],0,1,color="gray",label=newmark_label,linewidth=2,linestyle=:dash)
    Plots.plot!(sol_eigen(t)...,0,1,color="black",label=eigen_label,linewidth=2,linestyle=:dash)
    Plots.plot!(xlims=[2,6],ylims=[2,7])
    Plots.savefig("presentation/gifs/eigen"*@sprintf("/frame%i.png", j))
end

Plots.gif(anim, "img/framework/eigen.gif", fps=15)
