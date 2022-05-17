include("Beam1D.jl")
import Plots

x_grid = collect(0:0.1:1)
q(x,t) = -100*(0 ≤ x ≤ 0.6)*(0.4 ≤ x ≤ 1)
EI(x) = 1
mu(x) = 1

BC = [1,-1,0,0]
par = Beam1D.Parameters(mu,EI,q,BC)
sys = Beam1D.build(x_grid,par)

times = collect(0:0.05:1)
u = sys.S\Beam1D.evaluate(sys.f,0)
zero = zeros(length(u))
IC = [u zero zero]

sys = Beam1D.build(x_grid,par)
sol = Beam1D.solve_tr(sys, IC, times)

anim = Plots.@animate for j ∈ 1:length(times)
    Plots.plot(sol[j], LinRange(0:1e-3:1), 
        linewidth=:3, color="black")
end

Plots.gif(anim, "../Images/beam_animation.gif", fps=10)