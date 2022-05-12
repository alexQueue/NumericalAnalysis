include("Beam1D.jl")
import Plots
using DataStructures

x_0 = 1
xprime_0 = 0
M_0 = nothing
Q_0 = nothing

x_L = nothing
xprime_L = nothing
M_L = 0
Q_L = 0


# Parameter validations
nothing_count = counter([x_0, xprime_0, M_0, Q_0, x_L, xprime_L, M_L, Q_L])[nothing]
nothing_count == 4 ? nothing : throw(AssertionError("Must have exactly 4 BCs"))

BoundaryConditions = Beam1D.BoundaryConditions(
  x_0, xprime_0, M_0, Q_0,
  x_L, xprime_L, M_L, Q_L
)

x_grid = collect(0:0.01:1)
q = function(x)
    if 0.4 ≤ x ≤ 0.41
        -1
    elseif 0.89 ≤ x ≤ 0.9
        -1
    else
        0
    end
end

par = Beam1D.Parameters(x->1,x->1,q,BoundaryConditions)
sys = Beam1D.build(x_grid,par)
sol = Beam1D.solve_st(sys)

Plots.plot(sol)
Plots.plot!(sol,sys.x,seriestype=:scatter)
