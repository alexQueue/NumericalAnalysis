import Plots; Plots.pgfplotsx() 
using LaTeXStrings
include("Beam2D.jl")

function main(args)
    file = "data/"*args[1]
    problem = Beam2D.Problem(file*".ne")
    xy = Beam2D.u_to_Vh(problem, zeros(problem.size))
    p = Plots.plot(xy[1],xy[2],0,1,color="black",label=false,linewidth=2)
    Plots.plot!(p,title=L"\textcolor{blue}{FREE}\hspace{0.3cm}\textcolor{black}{FIXED} \hspace{0.3cm}\textcolor{red}{FORCE} \hspace{0.3cm}\textcolor{green}{MOVABLE}")

    colors = Dict("FREE" => "blue", "FIXED" => "black", "FORCE" => "red", "MOVABLE" => "green")
    for node in problem.nodes
        color = colors[node.type]
        Plots.scatter!(p,[node.coord[1]],[node.coord[2]],label=false,color=color,markersize=10)
    end
    Plots.savefig(p,file*".svg")
end

main(ARGS)
