include("Beam1D.jl")
import Plots
using DataStructures

L_0 = 0.0
L = 1.0
h = 0.01

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
BC = Dict("x_0"=>1,"xprime_0"=>0,"M_0"=>0,"Q_L"=>0)
# BC = Dict("x_L"=>1,"xprime_L"=>2,"M_L"=>3,"Q_L"=>4)
# BC = Dict("x_0"=>1,"xprime_0"=>2,"M_0"=>3,"Q_0"=>4)
BoundaryConditions = Beam1D.make_BC_from_dict(BC)

x_grid = collect(L_0:h:L)
q(x) = 5#-5*(L_0 ≤ x ≤ 0.6)*(0.4 ≤ x ≤ L)
E(x) = 1#x^2+1
I(x) = 1#x+1

par = Beam1D.Parameters(E,I,q,BoundaryConditions)
sys = Beam1D.build(x_grid,par)
sol = Beam1D.solve_st(sys)

u = sys.S\sys.f
println(u[1:4])


# Diagnostics
using LinearAlgebra
println("-------------")
println("At position 0:")
println(u[1:4])
println(string("X is ", u[1]))
println(string("X' is ", u[2]))
println(string("X'' is ", dot(u[1:4], -[6/h^2,   4/h,     -6/h^2,  2/h])))
println(string("X''' is ", dot(u[1:4], -[12/h^3,   6/h^2,   -12/h^3, 6/h^2])))
println("-------------")
println("At position L:")

println(string("X is ", u[end-1]))
println(string("X' is ", u[end]))
println(string("X'' is ", dot(u[end-3:end], [6/h^2,   2/h,  -6/h^2,     4/h])))
println(string("X''' is ", dot(u[end-3:end], [12/h^3,     6/h^2, -12/h^3,     6/h^2])))

# print(sol[1:4])

Plots.plot(sol)
Plots.plot!(sol,sys.x,seriestype=:scatter)