include("Beam2D.jl")
import Plots

problem = Beam2D.Problem("data/framework3.ne")
sys  = Beam2D.System(problem)

times = collect(LinRange(0,10,500))
Beam2D.vibrate_frame(sys, times, "img/framework3.gif")
