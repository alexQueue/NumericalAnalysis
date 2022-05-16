include("Beam1D.jl")
import Plots

x_grid = collect(0:0.01:1)
q(x) = -(0 ≤ x ≤ 0.6)*(0.4 ≤ x ≤ 1)

par = Beam1D.Parameters(x->1,x->1,q,[0,0,0,0])
sys = Beam1D.build(x_grid,par)
sol = Beam1D.solve_st(sys)

Plots.plot(sol)
Plots.plot!(sol,sys.x,seriestype=:scatter)