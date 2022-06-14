include("Beam1D.jl")
include("Diagnostics.jl")
include("analytic_testcases.jl")
import Plots
import LinearAlgebra

L = 1.0
q0 = 3
M_0 = 9
a = 1
EI_const = 1
#analytic solution with BCs and q function
ana_sol, BCs, q_func = case_7_constant_load(L , q0 ,EI_const)
""" just copy and paste to run example
ana_sol, BCs, q_func = case_1_supported_beam_(L , a, q0 ,EI_const)
#case_7_constant_load(L , q0 ,EI_const)
#case_8_partly_constant_load(L, q0, a,  EI_const)
#case_9_decreasing_load(L, q0,  EI_const)
#case_10_free_momentum_at_a(L, a, M_0, EI_const)
"""

#parameters 
pars = (mu=x->1,EI=x->EI_const,q=q_func)

resolutions = [3, 5, 11, 21, 51, 101, 201]
errors = []

for resolution in resolutions
  grid = collect(LinRange(0,L,resolution))
  sys  = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
  sol = Beam1D.solve_st(sys)
  error = LinearAlgebra.norm(ana_sol(grid)-sol(grid)) / resolution
  append!(errors, [error])
end
#plot
plt = Plots.plot(resolutions, errors, xaxis=:log, yaxis=:log, seriestype=:scatter, label = "Convergence")
Plots.ylabel!("error")
Plots.xlabel!("resolution")
println("Plotting convergence. Press enter to continue.")
Plots.gui(plt)
readline()

# NOTE: This should converge but it currently fails to.
