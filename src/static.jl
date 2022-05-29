include("Beam1D.jl")
include("Diagnostics.jl")
import Plots

L_0 = 0.0
L = 1.0
h = 0.01
# BC using dictionary instead
BC = Dict("x_0"=>1,"xprime_0"=>2,"M_0"=>3,"Q_0"=>4)
BoundaryConditions = Beam1D.make_BC_from_dict(BC)

x_grid = collect(L_0:h:L)
q(x,t) = (0 ≤ x ≤ 0.6)*(0.4 ≤ x ≤ 1)
EI(x) = x^2+1
μ(x) = x+1

par = Beam1D.Parameters(μ,EI,q,BoundaryConditions)

sys = Beam1D.build(x_grid,par)
sol, u = Beam1D.solve_st(sys)

#Diagnostics
Diagnostics.print_diagnostics(u, h)

plt = Plots.plot(sol)
Plots.plot!(sol,sys.x,seriestype=:scatter)
println("Plotting. Press enter to continue.")
Plots.gui(plt)
readline()