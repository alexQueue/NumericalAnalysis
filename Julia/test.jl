include("Beam1D.jl")
import Plots
using DataStructures

x_0 = 1
xprime_0 = 0
M_0 = nothing
Q_0 = nothing

x_L = 1
xprime_L = nothing
M_L = nothing
Q_L = 1


# Parameter validations
nothing_count = counter([x_0, xprime_0, M_0, Q_0, x_L, xprime_L, M_L, Q_L])[nothing]
nothing_count == 4 ? nothing : throw(AssertionError("Must have exactly 4 BCs"))

BoundaryConditions = Beam1D.BoundaryConditions(
  x_0, xprime_0, M_0, Q_0,
  x_L, xprime_L, M_L, Q_L
)

# BC using dictionary instead
BC = Dict("x_0"=>1,"xprime_0"=>0,"x_L"=>1,"Q_L"=>1)
BoundaryConditions = Beam1D.make_BC_from_dict(BC)

x_grid = collect(0:0.01:1)
q(x) = -(0 ≤ x ≤ 0.6)*(0.4 ≤ x ≤ 1)
E(x) = x^2+1
I(x) = x+1

par = Beam1D.Parameters(E,I,q,BoundaryConditions)
sys = Beam1D.build(x_grid,par)
sol = Beam1D.solve_st(sys)

Plots.plot(sol)
Plots.plot!(sol,sys.x,seriestype=:scatter)
