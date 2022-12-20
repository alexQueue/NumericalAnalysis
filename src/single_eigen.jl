include("Beam1D.jl")
import Plots
import LinearAlgebra
import LaTeXStrings

function get_eigen_cantilever(system::Beam1D.System, pars, L)
    n = 4
    A = L^(-1/2)
    x_j = [(x-0.5)*pi/L for x in 1:n]
    k = [x_j[x]/L for x in 1:n ]

    # analytic case: Assumption that EI and mu are constant
    freq = (pars.EI(0)./pars.mu(0)) .* k.^4
    w_j(x) =  [  
        A .* (
            (cosh.(k[j].*x) .- cos.(k[j] .*x)) 
            .- ((cosh.(x_j[j]) .- cos.(x_j[j]))./(sinh.(x_j[j]) .+ sin.(x_j[j])) 
            .* (sinh.(k[j] .*x) .- sin.(k[j] .*x)))
            ) for j in 1:n
        ]

    eigenvec = w_j(system.problem.grid)
    display(eigenvec)
    return freq, eigenvec
end

#TODO not complete
function solve_dy_eigen_analytical(system::Beam1D.System, pars, L)
    freqs, evecs = get_eigen_cantilever(system, pars, L )
    
    function get_sol(IC::Matrix{Float64})
        @assert size(IC) == (sys.shape[2],2) "Wrong IC size for given system"

        as = evecs\IC[:,1]
        bs = (evecs\IC[:,2])./freqs

        sol(t) = u_to_Vh(sys.problem.grid, evecs*(as.*cos.(freqs.*t).+bs.*sin.(freqs.*t)))

        return sol
    end
    return get_sol
end

#TODO: Check numerical freqs/modes against analytical for cantilever

#Initial condition
pars = (mu=x->1,EI=x->1,q=x->10)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'M')=>0,
            (1,'Q')=>0)
L = 1
grid = collect(LinRange(0,1,20))

prob = Beam1D.Problem(pars,BCs,grid)
sys  = Beam1D.System(prob)

IC   = [sys.Se\sys.qe zeros(sys.shape[2])]

#Change load & solve dynamically with eigen method
ts = collect(LinRange(0,10,100))

#get_sol = Beam1D.solve_dy_eigen(sys)
#sol     = get_sol(IC)

#anim = Plots.@animate for t in ts
#    Plots.plot([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
#    Plots.plot!(sol(t)...,0,1,ylim=[-2,2],color="black",label=false,linewidth=2)
#end

#Plots.gif(anim, "img/single/beam_animation_eigen_num.gif", fps=15)

#TODO compare numerical and analytical 

#Plot analytical solution 
freq, w_j = get_eigen_cantilever(sys, pars, L)

#plt = Plots.plot(w_j(grid), xlim = [0, 1])
#Plots.plot!([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
#Plots.savefig(plt, "img/single/eigenfunctions_analytical_test_pia.svg")

get_sol = solve_dy_eigen_analytical(sys, pars, L)
sol = get_sol(IC)

anim2 = Plots.@animate for t in ts
   Plots.plot([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
    Plots.plot!(sol(t)...,0,1,ylim=[-2,2],color="black",label=false,linewidth=2)
end

Plots.gif(anim2, "img/single/beam_animation_eigen_ana_PiaTest.gif", fps=15)



#numerical 
#evals, evecs, freqs, modes = Beam1D.get_vibrations(sys)

# Plot analytical against numerical eigenfrequencies/i.e. eigemvalues
#Plots.plot(freq[1:4], seriestype=:scatter, label = "analytical", mc="blue")
#Plots.plot!(freqs[1:4], seriestype=:scatter, label = "numerical", mc="red")
#Plots.xlabel!(LaTeXStrings.L"j")#
#Plots.ylabel!(LaTeXStrings.L" Eigenfrequency  $\omega_j$")
#Plots.plot!(legend=:topleft)
#Plots.savefig("img/single/Eigenfrequencies_analytical_numerical.svg")

#Generate some movies of eigemodes and superpoistions

#get_sol = Beam1D.solve_dy_eigen(sys)
#sol     = get_sol(IC)

#anim = Plots.@animate for t in ts
#    Plots.plot([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
#    Plots.plot!(sol(t)...,0,1,ylim=[-4,4],color="black",label=false,linewidth=2)
#end

#Plots.gif(anim, "img/single/beam_animation_eigen.gif", fps=15)



#TODO: Implement some other configurations for gifs
