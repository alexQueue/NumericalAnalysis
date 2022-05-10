include("Beam1D.jl")
import Plots
using DataStructures

# parameters:
E = 1
I = 1
q = -1
x_values = [0,0.2,1]


x_0 = 1
xprime_0 = 0
M_0 = nothing
Q_0 = nothing

x_L = nothing
xprime_L = nothing
M_L = 0
Q_L = 0

# TODO: Should be exactly 4 "nothings" in our BCs
nothing_count = counter([x_0, xprime_0, M_0, Q_0, x_L, xprime_L, M_L, Q_L])[nothing]
nothing_count == 4 ? nothing : throw(AssertionError("Must have exactly 4 BCs"))

BoundaryConditions = Beam1D.BoundaryConditions(
  x_0, xprime_0, M_0, Q_0,
  x_L, xprime_L, M_L, Q_L
)



par = Beam1D.Parameters(E, I, q, BoundaryConditions)
sys = Beam1D.build(x_values,par)

# sys_with_bcs = Beam1D.add_boundary_conditions(sys, BCs, par)

sol = Beam1D.solve_st(sys)

Plots.plot(sol)
Plots.plot!(sol,sys.x,seriestype=:scatter)


# testing...
# print(Beam1D.BCs_from_Moments(0, 0, par, sys))