include("Beam1D.jl")
import Plots
using LinearAlgebra

L_0 = 0.0
L = 1.0
h = 0.01
# BC using dictionary instead
BC = Dict("x_0"=>1,"xprime_0"=>2,"M_0"=>3,"Q_0"=>4)
BoundaryConditions = Beam1D.make_BC_from_dict(BC)

x_grid = collect(L_0:h:L)
q(x) = (0 ≤ x ≤ 0.6)*(0.4 ≤ x ≤ 1)
EI(x) = x^2+1
mu(x) = x+1

par = Beam1D.Parameters(mu,EI,q,BoundaryConditions)

sys = Beam1D.build(x_grid,par)
sol, u = Beam1D.solve_st(sys)

#Diagnostics
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


p = Plots.plot(sol)
Plots.plot!(sol,sys.x,seriestype=:scatter)

display(p)