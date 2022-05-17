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
		mu::Function
		EI::Function
		q::Function
		BCs::BoundaryConditions
	end

	struct System
		par::Parameters
		x::Vector{Float64}
		S::SparseArrays.SparseMatrixCSC{Float64,Int64}
		M::SparseArrays.SparseMatrixCSC{Float64,Int64}
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

	β = 1/4; γ = 1/2
	"""
	Newmark methdod for beams with initial conditions
	IC[1,2,3] = [u₀, u̇₀, ü₀] and M and S being
	the mass and stiffness matrices and p being the 
	forcing terms together with the boundary conditions.
	"""
	function solve_tr(
			sys::System, 
			IC::Matrix{Float64}, 
			times::Vector{Float64})

		nₓ = length(IC[:,1])
		nₜ = length(times)
		
		q(t) = zeros(nₓ)

		u = zeros(nₓ,nₜ); u̇ = zeros(nₓ,nₜ); ü	= zeros(nₓ,nₜ);
		u[:,1] = IC[:,1]
		u̇[:,1] = IC[:,2]
		ü[:,1] = IC[:,3]
		for j=1:nₜ-1
			hⱼ = times[j+1] - times[j]
			uⱼ_star = u[:,j] + u̇[:,j]*hⱼ + (1/2 - β)*ü[:,j]*hⱼ^2
			u̇ⱼ_star = u̇[:,j] + (1 - γ)*ü[:,j]*hⱼ
	
			Suⱼ_star = sys.S*uⱼ_star
			N_u = length(sys.x)*2
			i = [1,2,N_u,N_u-1]
			Suⱼ_star[i] .= 0

			ü[:,j+1] = (sys.M+β*hⱼ^2*sys.S)\(q(times[j+1]) - Suⱼ_star)
			u̇[:,j+1] = u̇ⱼ_star + γ*ü[:,j+1]*hⱼ
			u[:,j+1] = uⱼ_star + β*ü[:,j+1]*hⱼ^2
		end
		return [x -> CubicHermiteSpline.CubicHermiteSplineInterpolation(
			sys.x, u[1:2:end,j], u[2:2:end,j])(x) for j ∈ 1:length(times)]
	end

	# Applies Dirichlet and Neumann BCs at x=0 and x=L
	function set_x_BCs(S::SparseArrays.SparseMatrixCSC{Float64,Int64}, f::Vector{Float64}, N_u::Integer, par::Parameters)
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
	end

	function build(x::Vector{Float64},par::Parameters)
		#Shape functions and second derivatives on [0,h]
		phi_0(h,p) = [ (p/h)^2*(2*p/h-3)+1;
		               p*(p/h-1)^2        ;
		              -(p/h)^2*(2*p/h-3)  ;
		               p^2/h*(p/h-1)      ]
		phi_2(h,p) = [ (12*p-6*h)/h^3     ;
		               ( 6*p-4*h)/h^2     ;
		              -(12*p-6*h)/h^3     ;
		               ( 6*p-2*h)/h^2     ]

		#3 point Gaussian quadrature on [0,h]
		GQ3(h,f) = 5*h/18*f(h/2*(1-sqrt(3/5))) +
		           4*h/9 *f(h/2              ) +
		           5*h/18*f(h/2*(1+sqrt(3/5)))

		#Local System using 3 point Gaussian quadrature
		i_loc       = [1, 2, 3, 4]
		M_loc(h,p0) = GQ3(h,p->phi_0(h,p)*phi_0(h,p)'*par.mu(p+p0))
		S_loc(h,p0) = GQ3(h,p->phi_2(h,p)*phi_2(h,p)'*par.EI(p+p0))
		f_loc(h,p0) = GQ3(h,p->phi_0(h,p)*par.q(p+p0))
		
		#Global Variables
		N_v = length(x) #Number of vertices
		N_e = N_v-1     #Number of elements
		N_u = N_v*2     #Number of unknowns
		N_bc = 8		# Number of (possible) Boundary Conditions
		L = x[end]

		#Global System
		M = SparseArrays.spzeros(Float64,N_u,N_u)
		f = zeros(Float64,N_u + N_bc)
		S = SparseArrays.spzeros(Float64,N_u + N_bc,N_u)

		#Element contributions
		for k in 1:N_e
			i       = i_loc.+2*(k-1)
			h       = x[k+1]-x[k]
			
			M[i,i] += M_loc(h,x[k])
			S[i,i] += S_loc(h,x[k])
			f[i]   += f_loc(h,x[k])
		end

		#Boundary Conditions
		# <<<<<<< HEAD
		# =======
		# 		i      = [1,2,N_u-1,N_u] #Boundary indices
		# 		S[i,:] = SparseArrays.sparse(LinearAlgebra.I,N_u,N_u)[i,:]
		# 		M[i,:] .= 0
		# 		f[i]   = par.BCs
		# >>>>>>> main

		# S embeds boundary conditions in an 8xN block at the bottom of the matrix.
		# First four rows are for x and x' BCs.

		set_x_BCs(S, f, N_u, par)

		h_0 = x[2]-x[1]
		h_L = x[end]-x[end-1]

		if par.BCs.Q_0 !== nothing
			S[N_u + 5, 1:4] = Q_at_point(0.0, h_0, par)
			f[N_u + 5] = -par.BCs.Q_0 
		end
		if par.BCs.M_0 !== nothing
			S[N_u + 6, 1:4] = [6/h_0^2	 4/h_0	-6/h_0^2	 2/h_0] * 2/h_0 * par.EI(0)
			f[N_u + 6] = -par.BCs.M_0 
		end
		
		if par.BCs.Q_L !== nothing
			S[N_u + 7, end-3:end] = Q_at_point(L, h_L, par)
			f[N_u + 7] = -par.BCs.Q_L 
		end
		if par.BCs.M_L !== nothing
			S[N_u + 8, end-3:end] = [6/h_L^2	 2/h_L	-6/h_L^2	 4/h_L] * 2/h_L * par.EI(L)
			f[N_u + 8] = par.BCs.M_L 
		end
		
		#Packaging
		return System(par,x,S,M,f)
	end

	# Apply Q at a given point . Usually 0 or L for BCs
	# Returns a vector to place in our S matrix
	function Q_at_point(point::Float64, h::Float64, par::Parameters)		
		S_vector = [6/h^2,	 4/h,	-6/h^2,	 2/h] * 2/h * ForwardDiff.derivative(par.EI, point)
		S_vector += [12/h^3,	 6/h^2,	-12/h^3,	 6/h^2] * par.EI(point)
		
		return S_vector
	end
end