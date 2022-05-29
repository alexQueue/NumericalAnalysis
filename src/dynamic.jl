include("Beam1D.jl")
include("Diagnostics.jl")
import Plots

h = 0.05
x_grid = collect(0:h:1)
q(x,t) = 20*(0 ≤ x ≤ 0.6)*(0.4 ≤ x ≤ 1) + t*x
EI(x) = 1
mu(x) = 1


BC = Dict("x_0"=>1,"xprime_0"=>0,"M_L"=>0,"Q_L"=>0)
BC = Dict("x_0"=>1,"xprime_0"=>0,"xprime_L"=>0,"x_L"=>1)
BoundaryConditions = Beam1D.make_BC_from_dict(BC)

par = Beam1D.Parameters(mu,EI,q,BoundaryConditions)
sys = Beam1D.build(x_grid,par)

times = collect(0:0.4:2)
u = sys.S\Beam1D.evaluate(sys.f, 0)
zero = zeros(length(u))
IC = [u zero zero]

sys = Beam1D.build(x_grid,par)
sol, ut = Beam1D.solve_tr(sys, IC, times)

anim = Plots.@animate for j ∈ 1:length(times)
    Plots.plot(sol[j],LinRange(0:1e-3:1), linewidth=:3, color="black", ylim=(0, 2))
    Plots.plot!(sol[j],sys.x,seriestype=:scatter)
end

# TODO: Save to the same location regardless of where the file is run from.
Plots.gif(anim, "img/beam_animation.gif", fps=10)

Diagnostics.print_diagnostics(ut[:, end], h)
