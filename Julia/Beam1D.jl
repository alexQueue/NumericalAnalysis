module Beam1D
	import SparseArrays, LinearAlgebra, CubicHermiteSpline

	mutable struct Parameters
		E::Function
		I::Function
		q::Function
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

		# Non-constant load
		phi0(h) = [ 12/h^3	 8/h^2	-12/h^3	 4/h^2
					8/h^2	 16/3/h	-8/h^2	 8/3/h
				   -12/h^3	-8/h^2	 12/h^3	-4/h^2
					4/h^2	 8/3/h	-4/h^2	 4/3/h]

		phih2(h) = [0	 0		 0		 0
					0	 4/3/h	 0		-4/3/h
					0	 0		 0		 0
					0	-4/3/h	 0		 4/3/h]

		phih(h) = [ 12/h^3	 4/h^2	-12/h^3	 8/h^2
					4/h^2	 4/3/h	-4/h^2	 8/3/h
	   			    -12/h^3	-4/h^2	 12/h^3	-8/h^2
					8/h^2	 8/3/h	-8/h^2	 16/3/h]

		phi(h) = [h/3	 2h/3	 0
				  0		 h^2/6	 0
				  0		 2h/3	 h/3
				  0		-h^2/6	 0]

		# Simpsons rule
		S_loc(h,EI) = h/3 * (phi0(h)*EI[1] + phih2(h)*EI[2] + phih(h)*EI[3])
		f_loc(h,q) = h/3 * (phi(h)*q)

		#Global Variables
		N_v = length(x) # Number of vertices
		N_e = N_v-1     # Number of elements
		N_u = N_v*2     # Number of unknowns

		#Global System
		f = zeros(Float64,N_u)
		S = SparseArrays.spzeros(Float64,N_u,N_u)

		#Element contributions
		for k in 1:N_e
			h     = x[k+1]-x[k]
			i     = i_loc.+2*(k-1)
			
			x_loc = [x[k]; (x[k]+x[k+1])/2; x[k+1]]
			EI    = par.E.(x_loc).*par.I.(x_loc)
			q     = par.q.(x_loc)

			# Non-constant load
			S[i,i] += S_loc(h, EI)
			f[i]   += f_loc(h, q)
		end

		#Boundary Conditions
		i      = [1,2,N_u-1,N_u] #Boundary indices
		S[i,:] = SparseArrays.sparse(LinearAlgebra.I,N_u,N_u)[i,:]
		f[i]   = par.BCs

		#Packaging
		return System(par,x,S,f)
	end
end