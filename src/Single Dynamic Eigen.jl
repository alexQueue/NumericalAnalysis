include("Beam1D.jl")
import Plots

function get_eigen_cantilever(system::System, pars, L)
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
pars = (mu=x->1,EI=x->1,q=x->-10)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'M')=>0,
            (1,'Q')=>0)
grid = collect(LinRange(0,1,10))

prob = Beam1D.Problem(pars,BCs,grid)
sys  = Beam1D.System(prob)

IC   = [sys.Se\sys.qe zeros(sys.shape[2],2)]

#Change load & solve dynamically with eigen method
times = collect(LinRange(0,10,500))
modes, weights = Beam1D.solve_dy_eigen(sys,IC)

anim = Plots.@animate for t in times
    Plots.plot(x -> modes(x)'weights(t),ylim=[-1.5,1.5])
end

Plots.gif(anim, "beam_animation.gif", fps=15)

#TODO: Implement some other configurations for gifs