include("Beam1D.jl")
import Plots

#Cantilever
pars = (mu=x->1, EI=x->1, q=x->-10)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'M')=>0,
            (1,'Q')=>0)
grid = collect(LinRange(0,1,10))

sys = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
sol = Beam1D.solve_st(sys)

Plots.plot(sol)
Plots.plot!(sol,sys.problem.grid,seriestype=:scatter)

#TODO: make a plot showing different loads / B.C.s