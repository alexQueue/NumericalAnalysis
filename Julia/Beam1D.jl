module Beam1D
	import SparseArrays, LinearAlgebra, CubicHermiteSpline, ForwardDiff

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
		E::Function
		I::Function
		q::Function
		μ::Function
		BCs::BoundaryConditions
	end

	Base.@kwdef struct System
		par::Parameters
		x::Vector{Float64}
		S::SparseArrays.SparseMatrixCSC{Float64,Int64}
		M::SparseArrays.SparseMatrixCSC{Float64,Int64} = nothing
		f::Vector{Float64}
	end

	"""
	Makes an object of BoundaryConditions struct from a dictionary
	of the boundary values
	"""
	function make_BC_from_dict(bc_dict)
		BC_array = Array{Union{Float64, Nothing}}(undef, 8)
		fields = string.(fieldnames(BoundaryConditions))
		for (field,i) ∈ zip(fields, 1:length(fields))
			BC_array[i] = haskey(bc_dict, field) ? bc_dict[field] : nothing
		end
		BoundaryConditions(BC_array...)
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
		phi0_S(h) = [ 12/h^3	 8/h^2	-12/h^3	 4/h^2
					8/h^2	 16/3/h	-8/h^2	 8/3/h
				   -12/h^3	-8/h^2	 12/h^3	-4/h^2
					4/h^2	 8/3/h	-4/h^2	 4/3/h]

		phih2_S(h) = [0	 0		 0		 0
					0	 4/3/h	 0		-4/3/h
					0	 0		 0		 0
					0	-4/3/h	 0		 4/3/h]

		phih_S(h) = [ 12/h^3	 4/h^2	-12/h^3	 8/h^2
					4/h^2	 4/3/h	-4/h^2	 8/3/h
	   			    -12/h^3	-4/h^2	 12/h^3	-8/h^2
					8/h^2	 8/3/h	-8/h^2	 16/3/h]

		phih3_M(h) = [49	14h		140		-28h
					14h		4h^2	40h		-8h^2
					140		40h		400		-80h
					-28h	-8h^2	-80h	16h^2]/729

		phi2h3_M(h) = [49	14h		140		-28h
					14h		4h^2	40h		-8h^2
					140		40h		400		-80h
					-28h	-8h^2	-80h	16h^2]/729

		phi(h) = [h/3	 2h/3	 0
				  0		 h^2/6	 0
				  0		 2h/3	 h/3
				  0		-h^2/6	 0]

		# Simpsons rule
		S_loc(h,EI) = h/3 * (phi0_S(h)*EI[1] + phih2_S(h)*EI[2] + phih_S(h)*EI[3])
		f_loc(h,q) = h/3 * (phi(h)*q)
		M_loc(h,μ) = 9h/8 * (phih3_M(h/3)*μ[1] + phi2h3_M(2h/3)*μ[2]) # 3/8 rule

		#Global Variables
		N_v = length(x) # Number of vertices
		N_e = N_v-1     # Number of elements
		N_u = N_v*2     # Number of unknowns
		N_bc = 8		# Number of (possible) Boundary Conditions
		L = x[end]

		#Global System
		f = zeros(Float64,N_u + N_bc)
		S = SparseArrays.spzeros(Float64,N_u + N_bc,N_u)
		
		non_zero_μ = (par.μ.(x) == zeros(length(x))) ? false : true
		if non_zero_μ
			M = SparseArrays.spzeros(Float64,N_u + N_bc,N_u)
		end

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

			if non_zero_μ
				x_locM 	= [(2x[k]+x[k+1])/3; (x[k]+2x[k+1])/3]
				μ	  	= par.μ.(x_locM)
				M[i,i] += M_loc(h, μ)
			end
		end

		#Boundary Conditions

		# S embeds boundary conditions in an 8xN block at the bottom of the matrix.
		# First four rows are for x and x' BCs.

		h_0 = x[2]-x[1]
		h_L = x[end]-x[end-1]

		# TODO: Clean up this horrible mess.

		if par.BCs.x_0 !== nothing
			S[N_u + 1, 1		] = 1
			f[N_u + 1] = par.BCs.x_0
		end
		if par.BCs.xprime_0 !== nothing
			S[N_u + 2, 2		] = 1
			f[N_u + 2] = par.BCs.xprime_0
		end
		if par.BCs.x_L !== nothing
			S[N_u + 3, end-1] = 1
			f[N_u + 3] = par.BCs.x_L
		end
		if par.BCs.xprime_L !== nothing
			S[N_u + 4, end	] = 1
			f[N_u + 4] = par.BCs.xprime_L
		end

		if par.BCs.Q_0 !== nothing
			Iprime = ForwardDiff.derivative(par.I, 0)
			Eprime = ForwardDiff.derivative(par.E, 0)
			derivative_factor = Eprime * par.I(0) + par.E(0) * Iprime
			S[N_u + 5, 1:4] = [6/h_0^2	 4/h_0	-6/h_0^2	 2/h_0] * 2/h_0 * derivative_factor
			S[N_u + 5, 1:4] += [12/h_0^3	 6/h_0^2	-12/h_0^3	 6/h_0^2] * par.E(0) * par.I(0)
			f[N_u + 5] = -par.BCs.Q_0 
		end
		if par.BCs.M_0 !== nothing
			S[N_u + 6, 1:4] = [6/h_0^2	 4/h_0	-6/h_0^2	 2/h_0] * 2/h_0 * par.E(0) * par.I(0)
			f[N_u + 6] = -par.BCs.M_0 
		end
		
		if par.BCs.Q_L !== nothing
			Iprime = ForwardDiff.derivative(par.I, L)
			Eprime = ForwardDiff.derivative(par.E, L)
			derivative_factor = Eprime * par.I(L) + par.E(L) * Iprime
			S[N_u + 7, end-3:end] = [6/h_L^2,	 4/h_L,	-6/h_L^2,	 2/h_L] * 2/h_L * derivative_factor
			S[N_u + 7, end-3:end] += [12/h_L^3,	 6/h_L^2,	-12/h_L^3,	 6/h_L^2] * par.E(L) * par.I(L)
			f[N_u + 7] = -par.BCs.Q_L 
		end
		if par.BCs.M_L !== nothing
			S[N_u + 8, end-3:end] = [6/h_L^2	 2/h_L	-6/h_L^2	 4/h_L] * 2/h_L * par.E(L) * par.I(L)
			f[N_u + 8] = par.BCs.M_L 
		end
		
		#Packaging
		if non_zero_μ
			return System(par=par,x=x,S=S,f=f,M=M)
		else
			return System(par=par,x=x,S=S,f=f)
		end
	end
end