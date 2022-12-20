module Dynamic2D
    include("Beam2D.jl")
    using Printf
    import Plots

    function animate_framework(ne_file, gif_savefile, 
            presentation_savedir, extremas=[-Inf,-Inf,Inf,Inf])
        problem = Beam2D.Problem(ne_file)
        sys  = Beam2D.System(problem)

        times = collect(LinRange(0,100,100))

        u = sys.Se\sys.qe
        IC = [u zeros(size(u)...,2)]

        sys.qe[:] .= 0

        xy = Beam2D.solve_dy_Newmark(sys,IC,times)
        undeformed = Beam2D.u_to_Vh(sys.problem,zeros(size(u)...))

        anim = Plots.@animate for (j,t) in enumerate(times)
            Plots.plot(undeformed[1],undeformed[2],0,1,color="black",label=false,linewidth=2
                ,linestyle=:dot,xlims=extremas[[1,3]],ylims=extremas[[2,4]])
            Plots.plot!(xy[j][1],xy[j][2],0,1,color="black",label=false,linewidth=2)
            Plots.savefig(presentation_savedir*@sprintf("/frame%i.png", j))
        end
        fps=15
        Plots.gif(anim, gif_savefile, fps=fps)
    end
end
