include("Beam1D.jl")
import Plots
using LinearAlgebra

h = 0.05
x_grid = collect(0:h:1)
q(x) = 20*(0 ≤ x ≤ 0.6)*(0.4 ≤ x ≤ 1)
EI(x) = 1
mu(x) = 1


BC = Dict("x_0"=>1,"xprime_0"=>0,"M_L"=>0,"Q_L"=>0)
BoundaryConditions = Beam1D.make_BC_from_dict(BC)

par = Beam1D.Parameters(mu,EI,q,BoundaryConditions)
sys = Beam1D.build(x_grid,par)

times = collect(0:0.05:10)
u = sys.S\sys.f
zero = zeros(length(u))
IC = [u zero zero]

sys = Beam1D.build(x_grid,par)
sol, ut = Beam1D.solve_tr(sys, IC, times)

anim = Plots.@animate for j ∈ 1:length(times)
    Plots.plot(sol[j],LinRange(0:1e-3:1), linewidth=:3, color="black", ylim=(0, 2))
    Plots.plot!(sol[j],sys.x,seriestype=:scatter)
end

# TODO: Save to the same location regardless of where the file is run from.
Plots.gif(anim, "Images/beam_animation.gif", fps=10)

u = ut[:, end]

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

