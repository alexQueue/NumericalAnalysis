include("Beam1D.jl")
import Plots

#Cantilever
pars = (mu=x->1, EI=x->10000, q=x->-0.01)
BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'H')=>0,
            (1,'G')=>0)
grid = collect(LinRange(0,200,10))

sys   = Beam1D.System(Beam1D.Problem(pars,BCs,grid))
xs,ys = Beam1D.solve_st(sys)

#TODO: make a plot showing different loads / B.C.s

p = Plots.plot()
Plots.plot!.(xs,ys,0,1)
p