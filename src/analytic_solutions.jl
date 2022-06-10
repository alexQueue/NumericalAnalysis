include("Beam1D_2.jl")
include("Diagnostics.jl")
include("analytic_testcases.jl")
import Plots
import LinearAlgebra

L = 1.0
q0 = 3
M_0 = 9
a = 1
EI_const = 1
grid = collect(LinRange(0,L,20))

#analytic solution with BCs and q function
ana_sol, BCs, q_func = case_1_supported_beam_(L , a, q0 ,EI_const)
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
sys  = Beam1D.build(Beam1D.Problem(pars,BCs,grid))
sol= Beam1D.solve_st(sys)

println("L2 norm: ", LinearAlgebra.norm(ana_sol(grid)-sol(grid)))

#plot
plt = Plots.plot(sol,sys.problem.grid,seriestype=:scatter, label = "numerical solution")
Plots.plot!(ana_sol, label = "analytical solution")
Plots.ylabel!("w(x)")
Plots.xlabel!("x")
Plots.savefig("AnalyticalvsNumerical_momentum.pdf")
println("Plotting. Press enter to continue.")
Plots.gui(plt)
readline()