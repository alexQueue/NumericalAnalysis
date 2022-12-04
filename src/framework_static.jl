module Static2D
    include("Beam2D.jl")
    import Plots

    function solve_and_plot(ne_file, savefile)
        problem = Beam2D.Problem(ne_file)
        sys = Beam2D.System(problem)

        p = Beam2D.plot_static(sys)
        Plots.savefig(p, savefile)
    end
end
