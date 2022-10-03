# This file needs to exist for standard package management
module NumericalAnalysis
    export Beam1D, Beam2D

    include("Beam1D.jl")
    include("Beam2D.jl")
end
