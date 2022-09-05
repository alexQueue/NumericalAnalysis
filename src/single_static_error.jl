include("Beam1D.jl")
import Plots, LinearAlgebra

foeppl(xi, alpha, n)=  ifelse.(xi .> alpha, (xi .-alpha).^n, 0 )

function case_1_supported_beam_(l , F, a ,EI)
	alpha = a/l
	beta  = (l-a)/l
	Xi(x) = x./l

	BCs  = Dict((0,'H')=>0,
	            (0,'G')=>(((F*l^2)/6)*(beta-beta^3))/EI,
	            (1,'H')=>0,
	            (1,'G')=>-(((F*l^2)/6)*(alpha-alpha^3))/EI)

	q_func(x) = q0#ifelse.(x.== a, q0, 0) #TODO we need F and not q
	return  x-> ((((F*l^3)/6).*(beta*Xi(x) .* (1-beta^2 .- Xi(x).^2) .+ foeppl(Xi(x), alpha, 3)))/EI), BCs, q_func 
end

function case_7_constant_load(l , q0 ,EI)
	q_func(x) = q0
	Xi(x) = x./l 

	BCs  = Dict((0,'H')=>0,
	            (0,'G')=>0,
	            (1,'M')=>0,
	            (1,'Q')=>0)

	return  x-> (((q0*l^4)/24)* (6*Xi(x).^2 - 4*Xi(x).^3 + Xi(x).^4 ))/EI, BCs, q_func
end

function case_8_partly_constant_load(l, q0, a,  EI)
	q_func(x) = ifelse.(x.>= a, q0, 0)
	BCs  = Dict((0,'H')=>0,
	            (0,'G')=>0,
	            (1,'M')=>0,
	            (1,'Q')=>0)

	Xi(x) = x./l 
	alpha = a/l
	beta  = (l-a) /l
	return x-> ((q0*l^4)/24)* (foeppl(Xi(x), alpha, 4) - 4*beta*Xi(x).^3+ 6*beta*(2-beta)*Xi(x).^2)/EI, BCs, q_func
end

function case_9_decreasing_load(l , q0 ,EI)
	q_func(x) = q0-q0*x
	BCs  = Dict((0,'H')=>0,
	            (0,'G')=>0,
	            (1,'M')=>0,
	            (1,'Q')=>0)
	Xi(x) = x./l 
	return  x->  (q0*l^4)/120 * (10*Xi(x).^2-10*Xi(x).^3 + 5*Xi(x).^4 -Xi(x).^5)/EI, BCs, q_func
 end

function case_10_free_momentum_at_a(l, M_0, a, EI)
	q_func(x) = 0
	BCs  = Dict((0,'H')=>0,
	            (0,'G')=>0,
	            (1,'M')=>M_0,
	            (1,'Q')=>0)
	Xi(x) = x./l 
	alpha = a/l
	return x-> ((-M_0*l^2)/2 .* (Xi(x).^2- foeppl(Xi(x),alpha, 2)))/EI, BCs, q_func
end

L = 1.0
q0 = 3
M0 = 9
a = 1
EI_const = 1
grid = collect(LinRange(0,L,20))

# Analytic solution with BCs and q function
""" just copy and paste to run example
# analytic_sol, BCs, q_func = case_7_constant_load(L , q0 ,EI_const)
# case_1_supported_beam_(L , q0, a ,EI_const)
# case_7_constant_load(L , q0 ,EI_const)
# case_8_partly_constant_load(L, q0, a,  EI_const)
# case_9_decreasing_load(L, q0,  EI_const)
# case_10_free_momentum_at_a(L, M_0, a, EI_const)
"""

cases = Dict(
    "supported_beam"        => Dict("fnc"=>case_1_supported_beam_,      "input"=>(L,q0,a,EI_const)),
    "constant_load"         => Dict("fnc"=>case_7_constant_load,        "input"=>(L,q0,  EI_const)),
    "partly_constant_load"  => Dict("fnc"=>case_8_partly_constant_load, "input"=>(L,q0,a,EI_const)),
    "decreasing_load"       => Dict("fnc"=>case_9_decreasing_load,      "input"=>(L,q0,  EI_const)),
    "free_momentum_at_a"    => Dict("fnc"=>case_10_free_momentum_at_a,  "input"=>(L,M0,a,EI_const)),
)

for (key,value) in cases
    local analytic_sol, Bcs, q_func = value["fnc"](value["input"]...)
    # Parameters 
    local pars = (mu=x->1 ,EI=x->EI_const, q=q_func)

    # Build and solve
    local problem = Beam1D.Problem(pars, BCs, grid)
    local sys = Beam1D.System(problem)


    local u_numeric = sys.Se\sys.qe # 4 boundary conditions at end
    local u_analytic = analytic_sol(grid)

    word_print = uppercasefirst(replace(key,"_"=>" "))
    println(word_print * ", L2 norm: \n\t", LinearAlgebra.norm(u_numeric[1:2:end-4]- u_analytic))

    local xs,ys = Beam1D.u_to_Vh(grid,u_numeric)

    # Plot
    local p = Plots.plot()
    Plots.plot!(grid,u_analytic, label = "Analytical solution",color="blue",legend=:topleft)
    Plots.plot!(Beam1D.eval(xs,0.5),Beam1D.eval(ys,0.5),seriestype=:scatter,label="Numerical Solution",color="red")
    Plots.ylabel!("w(x)")
    Plots.xlabel!("x")
    Plots.savefig("img/"*key*".png")
end

#TODO: Plot all the required configurations
