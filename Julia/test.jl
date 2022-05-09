include("Beam1D.jl")
import Plots

par = Beam1D.Parameters(1,1,-1,[0,0,0,0])
sys = Beam1D.build([0,0.2,1],par)
sol = Beam1D.solve_st(sys)

Plots.plot(sol)
Plots.plot!(sol,sys.x,seriestype=:scatter)