include("Beam2D.jl")
import Plots

param_texts = Dict(
    "E"=>"Young's modulus",
    "I"=>"Area moment of inertia",
    "A"=>"Stiffness",
)

function solve_with_param_change(param::String, only_end_beam::Bool=true)
    problem = Beam2D.Problem("data/crane.ne")

    funcs = [x->t for t in collect(LinRange(4,20,5))]

    p = Plots.plot()
    colors = ["red","green","blue","black","orange"]
    for (func,color) in zip(funcs,Iterators.cycle(colors))
        symbol = Symbol(param)
        start = ifelse(only_end_beam, length(problem.edges), 1)
        for i in start:length(problem.edges)
            setfield!(problem.edges[i], symbol, func)
        end

        sys = Beam2D.System(problem)
        u = sys.Se\sys.qe
        xs,ys = Beam2D.u_to_Vh(problem, u)

        label = "$param = $(func(0))"
        if only_end_beam
            Plots.plot!(p,xs[end],ys[end],0,1,label=label,
                xlims=(4.7,5.1), linewidth=2, color=color,
                legend=:topleft)
            Plots.plot!(p,title="$(param_texts[param]) effect on last beam")
        else
            Plots.plot!(p,xs[1],ys[1],0,1,label=label,color=color,linewidth=2)
            Plots.plot!(p,xs[2:end],ys[2:end],0,1,label=false,
                linewidth=2,color=color, xlims=(1.8,6),
                legend=:bottomright)
            Plots.plot!(p,title="$(param_texts[param]) effect on whole framework")
        end
    end
    if only_end_beam
        Plots.savefig(p, "img/framework/$(param)_comparison_end.svg")
    else
        Plots.savefig(p, "img/framework/$(param)_comparison.svg")
    end
end
