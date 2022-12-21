include("Beam1D.jl")
import Plots
import LinearAlgebra
import LaTeXStrings

function get_eigen_cantilever(system::Beam1D.System, pars, L, n=4)
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
    eigenvec = mapreduce(permutedims, vcat, eigenvec)
    return freq, eigenvec
end

#TODO not complete
function solve_dy_eigen_analytical(system::Beam1D.System, pars, L, n_m)
    freqs, evecs = get_eigen_cantilever(system, pars, L , n_m)
    
    function get_sol(IC::Matrix{Float64})
        @assert size(IC) == (sys.shape[2],2) "Wrong IC size for given system"
        
        midpoint_idx = Int32(round(size(evecs)[2]/2))
        #as = (L/6).* ((evecs[:, 1].*IC[1, 1])+ 4*(evecs[:, midpoint_idx].*IC[midpoint_idx, 1])+ (evecs[:, end].*IC[end-5, 1]))
        #bs = ((L/6).* ((evecs[:, 1].*IC[2, 1])+ 4*(evecs[:, midpoint_idx].*IC[midpoint_idx+1, 1])+ (evecs[:, end].*IC[end-4, 1])))./freqs_num
        

        as = evecs'\IC[1:2:end-4,1]
        bs = (evecs'\IC[2:2:end-4,2])./freqs
        #as_num = evecs_num[1:2:end-4,1]\IC[1:2:end-4,1]
        #bs_num = (evecs_num[1:2:end-4,1]\IC[2:2:end-4,2])./freqs_num

        sol(t) = evecs'*(as.*cos.(freqs.*t).+bs.*sin.(freqs.*t)) 
        #sol(t) = evecs_num[1:2:end-4, :]*(as_num.*cos.(freqs_num.*t).+bs_num.*sin.(freqs_num.*t)) #evecs'*(as.*cos.(freqs.*t).+bs.*sin.(freqs.*t)) 

        return sol
    end
    return get_sol
end

#Initial condition
pars = (mu=x->1,EI=x->1,q=x->10)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'Q')=>0,
            (1,'M')=>0)
L = 1
grid = collect(LinRange(0,1,200))

prob = Beam1D.Problem(pars,BCs,grid)
sys  = Beam1D.System(prob)

IC   = [sys.Se\sys.qe zeros(sys.shape[2])]

#Change load & solve dynamically with eigen method
ts = collect(LinRange(0,10,100))

get_sol_num = Beam1D.solve_dy_eigen(sys)
sol_num     = get_sol_num(IC)
get_sol_ana = solve_dy_eigen_analytical(sys, pars, L, 4)
sol_ana = get_sol_ana(IC)


anim = Plots.@animate for (j, t) in enumerate(ts)
    Plots.plot([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
    Plots.plot!(sol_num(t)...,0,1,ylim=[-2,2],color="black",label=false,linewidth=2)
    #Plots.plot!(grid,sol_ana(t),ylim=[-2,2],color="grey",label=false,linewidth=2)
    Plots.savefig("presentation/gifs/superpos"*@sprintf("/frame%i.png", j))
end

Plots.gif(anim, "img/single/beam_animation_eigen_num_ana.gif", fps=15)


#Plot analytical solution 
freq, w_j = get_eigen_cantilever(sys, pars, L)

#plt = Plots.plot(w_j(grid), xlim = [0, 1])
#Plots.plot!([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
#Plots.savefig(plt, "img/single/eigenfunctions_analytical_test_pia.svg")

get_sol = solve_dy_eigen_analytical(sys, pars, L)
sol = get_sol(IC)

anim2 = Plots.@animate for t in ts
  Plots.plot([0,1],[0,0],color="black",label=false,linewidth=2,linestyle=:dot)
  Plots.plot!(grid,sol(t),ylim=[-2,2],color="black",label=false,linewidth=2)
end

Plots.gif(anim2, "img/single/beam_animation_eigen_ana_PiaTest.gif", fps=15)


#numerical 
evals, evecs, freqs, modes = Beam1D.get_vibrations(sys)

# Plot analytical against numerical eigenfrequencies/i.e. eigemvalues
Plots.plot(freq[1:4], seriestype=:scatter, label = "analytical", mc="blue")
Plots.plot!(freqs[1:4], seriestype=:scatter, label = "numerical", mc="red")
Plots.xlabel!(LaTeXStrings.L"j")#
Plots.ylabel!(LaTeXStrings.L" Eigenfrequency  $\omega_j$")
Plots.plot!(legend=:topleft)
Plots.savefig("img/single/Eigenfrequencies_analytical_numerical.svg")

