include("Beam1D.jl")
import Plots

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

par = Beam1D.Parameters(x->1,x->1,q,[0,0,0,0])
sys = Beam1D.build(x_grid,par)
sol = Beam1D.solve_st(sys)

Plots.plot(sol)
Plots.plot!(sol,sys.x,seriestype=:scatter)