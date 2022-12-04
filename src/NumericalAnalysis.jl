# This file needs to exist for standard package management
module NumericalAnalysis
    export Beam1D, Beam2D
    export BasisFunctions
    export Dynamic2D, Static2D

    include("Beam1D.jl")
    include("Beam2D.jl")
    include("framework_dynamic.jl")
    include("framework_static.jl")
    include("basis_functions.jl")
end
