include("Beam1D.jl")
import Plots, LinearAlgebra

foeppl(xi, alpha, n)=  ifelse.(xi .> alpha, (xi .-alpha).^n, 0 )

function case_1_supported_beam_(l , a, F ,EI)
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

function case_10_free_momentum_at_a(l, a, M_0, EI)
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
M_0 = 9
a = 1
EI_const = 1
grid = collect(LinRange(0,L,20))

#analytic solution with BCs and q function
ana_sol, BCs, q_func = case_7_constant_load(L , q0 ,EI_const)
""" just copy and paste to run example
#case_1_supported_beam_(L , a, q0 ,EI_const)
#case_7_constant_load(L , q0 ,EI_const)
#case_8_partly_constant_load(L, q0, a,  EI_const)
#case_9_decreasing_load(L, q0,  EI_const)
#case_10_free_momentum_at_a(L, a, M_0, EI_const)

"""

#parameters 
pars = (mu=x->1,EI=x->EI_const,q=q_func)

# build and solve
sys  = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
sol= Beam1D.solve_st(sys)

println("L2 norm: ", LinearAlgebra.norm(ana_sol(grid)-sol(grid)))

#plot
plt = Plots.plot(sol,sys.problem.grid,seriestype=:scatter, label = "numerical solution")
Plots.plot!(ana_sol, label = "analytical solution")
Plots.ylabel!("w(x)")
Plots.xlabel!("x")
#Plots.savefig("AnalyticalvsNumerical_momentum.pdf")
println("Plotting. Press enter to continue.")
Plots.gui(plt)
readline()

#TODO: Plot all the required configurations