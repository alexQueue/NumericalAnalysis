include("Beam1D.jl")
import Plots
import LinearAlgebra

function get_eigen_cantilever(system::Beam1D.System, pars, L)
    n = 4
    A = L^(-1/2)
    x_j = [(x-0.5)*pi/L for x in 1:n]
    k = [x_j[x]/L for x in 1:n ]
    # analytic case: Assumption that EI and mu are constant
    freq = (pars.EI(0)./pars.mu(0)) .* k.^4
    w_j(x) =  [  
        x -> A .* (
            (cosh(k[j].*x) .- cos(k[j] .*x)) 
            .- ((cosh(x_j[j]) .- cos(x_j[j]))./(sinh(x_j[j]) .+ sin(x_j[j])) 
            .* (sinh(k[j] .*x) .- sin(k[j] .*x)))
            ) for j in 1:n
        ]

    return freq, w_j
end

#TODO: Check numerical freqs/modes against analytical for cantilever

#Initial condition
pars = (mu=x->1,EI=x->1,q=x->-100)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'H')=>0,
            (1,'G')=>0)
grid = collect(LinRange(0,1,20))

prob = Beam1D.Problem(pars,BCs,grid)
sys  = Beam1D.System(prob)

IC   = [sys.Se\sys.qe zeros(sys.shape[2])]

#Change load & solve dynamically with eigen method
ts = collect(LinRange(0,1,20))

get_sol = Beam1D.solve_dy_eigen(sys)
sol     = get_sol(IC)

anim = Plots.@animate for t in ts
    Plots.plot([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
    Plots.plot!(sol(t)...,0,1,ylim=[-1.5,1.5],color="black",label=false,linewidth=2)
end

Plots.gif(anim, "img/beam_animation_eigen.gif", fps=15)

#TODO: Implement some other configurations for gifs
