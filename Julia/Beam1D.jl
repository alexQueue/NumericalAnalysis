module Beam1D
	import SparseArrays, LinearAlgebra, CubicHermiteSpline

	mutable struct Parameters
		E::Float64
		I::Float64
		q::Float64
		BCs::Vector{Float64}
	end

	struct System
		par::Parameters
		x::Vector{Float64}
		S::SparseArrays.SparseMatrixCSC{Float64,Int64}
		f::Vector{Float64}
	end

	function solve_st(sys::System) #Stationary solver
		u = sys.S\sys.f
		return x -> CubicHermiteSpline.CubicHermiteSplineInterpolation(
		    sys.x,u[1:2:end],u[2:2:end])(x)
	end

	function solve_tr(sys::System) #Transient solver
		return 0
	end

	function build(x::Vector{Float64},par::Parameters)
		#Local System
		i_loc    = [ 1,     2,     3,     4  ]
		S_loc(h) = [ 6/h/h  3/h   -6/h/h  3/h;
		             3/h    2     -3/h    1  ;
		            -6/h/h -3/h    6/h/h -3/h;
		             3/h    1     -3/h    2  ]*2/h*par.E*par.I
		f_loc(h) = [ 1,     h/6,   1,    -h/6]*h/2*par.q

		#Global Variables
		N_v = length(x) #Number of vertices
		N_e = N_v-1     #Number of elements
		N_u = N_v*2     #Number of unknowns

		#Global System
		f = zeros(Float64,N_u)
		S = SparseArrays.spzeros(Float64,N_u,N_u)

		#Element contributions
		for k in 1:N_e
			h       = x[k+1]-x[k]
			i       = i_loc.+2*(k-1)
			S[i,i] += S_loc(h)
			f[i]   += f_loc(h)
		end

		#Boundary Conditions
		i      = [1,2,N_u-1,N_u] #Boundary indices
		S[i,:] = SparseArrays.sparse(LinearAlgebra.I,N_u,N_u)[i,:]
		f[i]   = par.BCs

		#Packaging
		return System(par,x,S,f)
	end
end