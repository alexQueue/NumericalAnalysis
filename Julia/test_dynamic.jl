include("Beam1D.jl")
import Plots

x_grid = collect(0:0.1:1)
t = 0:0.05:10
q(x) = -(0 ≤ x ≤ 0.6)*(0.4 ≤ x ≤ 1)
EI(x) = 1
mu(x) = 1

u = zeros(length(x_grid)*2)
n_unknowns = length(u)
IC = [u u u]

par = Beam1D.Parameters(mu,EI,q,[0,0,0,0],IC,t)
sys = Beam1D.build(x_grid,par)
sol = Beam1D.solve_tr(sys)

anim = Plots.@animate for j ∈ 2:length(sys.par.t)
    Plots.plot(sol(x_grid,j), ylim=(-2,2))
end

Plots.gif(anim, "beam_animation.gif")