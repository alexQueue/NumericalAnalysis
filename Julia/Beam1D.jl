module Beam1D
	import SparseArrays, LinearAlgebra, CubicHermiteSpline

	mutable struct BoundaryConditions
		x_0::Union{Float64, Nothing}
		xprime_0::Union{Float64, Nothing}
		M_0::Union{Float64, Nothing}
		Q_0::Union{Float64, Nothing}

		x_L::Union{Float64, Nothing}
		xprime_L::Union{Float64, Nothing}
		M_L::Union{Float64, Nothing}
		Q_L::Union{Float64, Nothing}
	end

	mutable struct Parameters
		E::Float64
		I::Float64
		q::Float64
		BCs::BoundaryConditions
	end

	struct System
		par::Parameters
		x::Vector{Float64}
		S::SparseArrays.SparseMatrixCSC{Float64,Int64}
		f::Vector{Float64}
	end

	# Simply solve u_1, u_2 (u3, u4?) from a given M_0
	function BCs_from_Moments(M_0, M_L, par::Parameters, sys::System)
		S_loc(h) = [ 6/h/h  3/h   -6/h/h  3/h;
		             3/h    2     -3/h    1  ;
		            -6/h/h -3/h    6/h/h -3/h;
		             3/h    1     -3/h    2  ]*2/h

		# TODO: E, I evaluated at a point.
		f = sys.f
		f[2] = -M_0 / par.E / par.I
		# f[1] = 
		f[end] = M_L / par.E / par.I
		# f[end - 1] = 

		u = sys.S\sys.f

		return u
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
		f = zeros(Float64,N_u + 8)
		S = SparseArrays.spzeros(Float64,N_u + 8,N_u)

		#Element contributions
		for k in 1:N_e
			h       = x[k+1]-x[k]
			i       = i_loc.+2*(k-1)
			S[i,i] += S_loc(h)
			f[i]   += f_loc(h)
		end

		#Boundary Conditions

		# S embeds boundary conditions in an 8xN block at the bottom of the matrix.
		# First four rows are for x and x' BCs.
		
		
		
		


		h_0 = x[2]-x[1]
		h_L = x[end]-x[end-1]

		# TODO: Clean up this horrible mess.
		if par.BCs.x_0 != nothing
			S[N_u + 1, 1		] = 1
			f[N_u + 1] = par.BCs.x_0
		end
		if par.BCs.xprime_0 != nothing
			S[N_u + 2, 2		] = 1
			f[N_u + 2] = par.BCs.xprime_0
		end
		if par.BCs.x_L != nothing
			S[N_u + 3, end-1] = 1
			f[N_u + 3] = par.BCs.x_L
		end
		if par.BCs.xprime_L != nothing
			S[N_u + 4, end	] = 1
			f[N_u + 4] = par.BCs.xprime_L
		end

		if par.BCs.Q_0 != nothing
			S[N_u + 5, 1:4] = [6/h_0^2	 4/h_0	-6/h_0^2	 2/h_0] * 2/h_0 * 0 # TODO: Not 0 but (E' * I + E * I')
			S[N_u + 5, 1:4] += [12/h_0^3	 6/h_0^2	-12/h_0^3	 6/h_0^2] * par.E * par.I
			f[N_u + 5] = -par.BCs.Q_0 
		end
		if par.BCs.M_0 != nothing
			S[N_u + 6, 1:4] = [6/h_0^2	 4/h_0	-6/h_0^2	 2/h_0] * 2/h_0 * par.E * par.I
			f[N_u + 6] = -par.BCs.M_0 
		end
		
		if par.BCs.Q_L != nothing
			S[N_u + 7, end-3:end] = [6/h_L^2,	 4/h_L,	-6/h_L^2,	 2/h_L] * 2/h_L * 0 # TODO: Not 0 but (E' * I + E * I')
			S[N_u + 7, end-3:end] += [12/h_L^3,	 6/h_L^2,	-12/h_L^3,	 6/h_L^2] * par.E * par.I
			f[N_u + 7] = -par.BCs.Q_L 
		end
		if par.BCs.M_L != nothing
			S[N_u + 8, end-3:end] = [6/h_L^2	 2/h_L	-6/h_L^2	 4/h_L] * 2/h_L * par.E * par.I
			f[N_u + 8] = par.BCs.M_L 
		end
		

		#Packaging
		return System(par,x,S,f)
	end
end