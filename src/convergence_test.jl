include("Beam1D.jl")
include("Diagnostics.jl")
include("analytic_testcases.jl")
import Plots
import LinearAlgebra

pars = (mu=x->1, EI=x->1, q=x->-(x<0.5))

BCs  = Dict((0,'H')=>1,
            (0,'G')=>2,
            (1,'M')=>3,
            (1,'Q')=>4)


resolutions = [10, 20, 50, 100, 200, 400, 800]

best_resolution = last(resolutions)*10
grid = collect(LinRange(0, 1, best_resolution))
sys   = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
best = Beam1D.solve_st(sys)


errors = []
for resolution in resolutions
  grid = collect(LinRange(0,1,resolution))
  sys  = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
  sol = Beam1D.solve_st(sys)
  error = LinearAlgebra.norm(best(grid) - sol(grid)) / resolution
  append!(errors, [error])
end
#plot
plt = Plots.plot(resolutions, errors, xaxis=:log, yaxis=:log, seriestype=:scatter, label = "Convergence", legend=:bottomleft)
Plots.ylabel!("error")
Plots.xlabel!("resolution")

ref1 = errors[1]*resolutions[1]*resolutions.^-1
ref1 = reshape(ref1, length(errors), 1)
ref2 = errors[1]*resolutions[1]^2*resolutions.^-2
ref2 = reshape(ref2, length(errors), 1)

Plots.plot!(resolutions, ref1, xaxis=:log, yaxis=:log, label = "Reference 1st Order")
Plots.plot!(resolutions, ref2, xaxis=:log, yaxis=:log, label = "Reference 2nd Order")

println("Plotting convergence. Press enter to continue.")

Plots.gui(plt)
readline()
