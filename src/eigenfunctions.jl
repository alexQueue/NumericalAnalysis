include("Beam1D.jl")
include("Beam2D.jl")
import CubicHermiteSpline
import LinearAlgebra
import Plots
using LaTeXStrings

function plot_eigenfunctions()
    #Initial condition
    pars = (mu=x->1,EI=x->1,q=x->-1)
    BCs  = Dict((0,'H')=>0,
                (0,'G')=>0,
                (1,'M')=>0,
                (1,'Q')=>0)
    grid = collect(LinRange(0,1,10))
    grid = collect(LinRange(0,1,40))

    #Change load & solve dynamically
    pars  = (mu=x->1, EI=x->1, q=x->-(x<0.5))
    sys  = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
    alpha = 2*[1, 1, 1, 1, 1, 1]#[4 , 20, 10, 50, 10]
    beta = alpha

    #compute eigenfunction and frequency numerically and analytically
    eigenvals, eigenvec, w = Beam1D.get_eigenvalues_numerical(sys, 4)
    eigenfreq, eigenfunc_ana = Beam1D.get_eigenvalues_analytical(sys, pars, 1, 4)

    #Plot eigenfrequencies (eigenvalues)
    plt = Plots.plot(eigenvals, seriestype=:scatter, label = "Numerical" )
    Plots.plot!(plt, eigenfreq, seriestype =:scatter, label = "Analytical")
    Plots.plot!(plt, plot_title="Eigenfrequencies", legend=:topleft)
    Plots.savefig(plt, "img/single/eigenfrequencies_analytical_numerical.svg")

    #Plot analytical eigenfunctions
    p1 = Plots.plot((eigenfunc_ana(grid))[1], xlim = [0, 1], title = L"w_1", color="blue", linewidth=2)
    p2 = Plots.plot((eigenfunc_ana(grid))[2], xlim = [0, 1], title = L"w_2", color="blue", linewidth=2)
    p3 = Plots.plot((eigenfunc_ana(grid))[3], xlim = [0, 1], title = L"w_3", color="blue", linewidth=2)
    p4 = Plots.plot((eigenfunc_ana(grid))[4], xlim = [0, 1], title = L"w_4", color="blue", linewidth=2)
    plt = Plots.plot(p1, p2, p3, p4, layout= (2, 2), legend= false)
    Plots.savefig(plt, "img/single/eigenfunctions_analytical.svg")

    #Plot numerical eigenfunctions
    pp1 = Plots.plot(
        Beam1D.u_to_Vh(grid, eigenvec[:, 1])..., 0, 1, title = L"w_1", yticks=LinRange(0.000,0.100,3), color="blue", linewidth=2)
    pp2 = Plots.plot(
        Beam1D.u_to_Vh(grid, eigenvec[:, 2])..., 0, 1, title = L"w_2", yticks=LinRange(-0.008,0.0045,3), color="blue", linewidth=2)
    pp3 = Plots.plot(
        Beam1D.u_to_Vh(grid, eigenvec[:, 3])..., 0, 1, title = L"w_3", yticks=LinRange(-0.001,0.002,3), color="blue", linewidth=2)
    pp4 = Plots.plot(
        Beam1D.u_to_Vh(grid, eigenvec[:, 4])..., 0, 1, title = L"w_4", yticks=LinRange(-0.0006,0.0006,3), color="blue", linewidth=2)
    plt = Plots.plot(pp1, pp2, pp3, pp4, layout= (2, 2), legend= false)
    Plots.savefig(plt, "img/single/eigenfunctions_numerical.svg")

    plt = Plots.plot(p1, pp1, p2, pp2, p3, pp3,  p4, pp4, layout=(4,2), legend=false)
    Plots.plot!(plt, plot_title="Analytical vs. Numerical eigenfunctions")
    Plots.savefig(plt, "img/single/eigenfunctions_numerical_and_analytical.svg")
   
    n = 4
    times = collect(LinRange(0,0.1,100))
    eigenfunc_t = [
        (x,t) -> CubicHermiteSpline.CubicHermiteSplineInterpolation(
        	sys.problem.grid,
            eigenvec[1:2:end-2, i],
            eigenvec[2:2:end-2, i])(x) .* real(w(t, alpha[i], beta[i])) for i in 1:n]

    numerical_anim = Plots.Animation()
    analytical_anim = Plots.Animation()
    for t in times
        res = zeros(length(grid))
        for i in 1:n
            res = res .+ eigenfunc_t[i](grid, t)
        end
        plt = Plots.plot(grid, res, ylim=[-1.5,1.5], xlim =[0,1],
            title="Superposition of eigenvalues, numerical",
            label = "w(x)", linewidth=2, color="blue")
        Plots.frame(numerical_anim, plt)
    end
    Plots.gif(numerical_anim, "img/single/beam_animation_eigenvals_superpos_numerical.gif", fps=15)
    
    anim = Plots.@animate for (j,t) in enumerate(times)
         #superposition
         res = zeros(length(grid))
         for i in 1:n
              res = res .+ eigenfunc_t[i](grid, t)
         end
         Plots.plot(grid, res, ylim=[-1.5,1.5], xlim =[0,1],
            title="Superposition of eigenvalues, numerical",
            label = "w(x)", linewidth=2, color="blue")
    end

    Plots.gif(anim, "img/single/beam_animation_eigenvals_superposition.gif", fps=15)  

    for i in 1:n 
        min_y = 1.5
        max_y = -1.5
        for t in times[1:end]
            vals = eigenfunc_t[i](grid, t)
            (min_v, max_v) = extrema(vals)
            if min_v < min_y
                min_y = min_v
            elseif max_v > max_y
                max_y = max_v
            end
        end
        anim = Plots.@animate for (j,t) in enumerate(times)
            Plots.plot(grid,
                title="Eigenvalue "*string(i),
                eigenfunc_t[i](grid, t),
                ylim=[min_y-0.3*abs(min_y),max_y+0.3*abs(max_y)], xlim =[0,1],
                linewidth=2, color="blue", label="w_"*string(i)*"(x)")
            end
        Plots.gif(anim,
            string("img/single/beam_animation_eigenvals_", i, ".gif"),
            fps=15
        ) 
    end
end
