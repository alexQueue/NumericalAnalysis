include("Beam2D.jl")
import Plots

function change_edge_param(edge::Beam2D.Edge, param::Function, param_name::String)
    func = Symbol(param_name)
    getfield(edge, func) = param
end

function solve_with_param_change()
    # problem = Beam2D.Problem("data/crane.ne")
    Es = [x->t for t in collect(LinRange(2.5,20,5))]
    p = Plots.plot()
    problem = Beam2D.Problem("data/crane.ne")
    colors = ["red","green","blue","black","orange"]
    for (E,color) in zip(Es,Iterators.cycle(colors))
        problem.edges[end].E = E

        sys = Beam2D.System(problem)
        u = sys.Se\sys.qe
        xs,ys = Beam2D.u_to_Vh(problem, u)
        Plots.plot!(p,xs,ys,0,1,label=string(E(0)),linewidth=2,legend=false,color=color)
    end
    Plots.savefig(p, "img/framework_E_comparison.svg")
end
