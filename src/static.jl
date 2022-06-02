include("Beam1D.jl")
include("Diagnostics.jl")
import Plots

pars = (mu=x->1, EI=x->1, q=x->-(x<0.5))

BCs  = Dict((0,'H')=>1,
            (0,'G')=>2,
            (1,'M')=>3,
            (1,'Q')=>4)
n = 10
grid = collect(LinRange(0,1,n))
h = grid[2] - grid[1]

sys   = Beam1D.build(Beam1D.Problem(pars,BCs,grid))
static_sol, u = Beam1D.solve_st(sys)

Diagnostics.print_diagnostics(u, h)

plt = Plots.plot!(static_sol, sys.problem.grid, seriestype=:scatter)
Plots.gui(plt)
println("Plotting. Press enter to continue.")
readline()

