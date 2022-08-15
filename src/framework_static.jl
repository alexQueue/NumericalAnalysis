include("Beam2D.jl")
import Plots

problem = Beam2D.Problem("data/crane.ne")
sys = Beam2D.System(problem)

p = Beam2D.plot_static(sys)
Plots.savefig(p, "img/framework_static.svg")
