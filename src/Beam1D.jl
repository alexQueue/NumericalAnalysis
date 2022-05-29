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

		indices::Array{Int64}
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
		f::Vector{Function}
	end

	"""
	Makes an object of BoundaryConditions struct from a dictionary
	of the boundary values
	"""
	function make_BC_from_dict(bc_dict)
		length(bc_dict) == 4 ? nothing : throw(AssertionError("Must have exactly 4 BCs"))
		BC_array = Array{Union{Float64, Nothing}}(undef, 8)
		fields = string.(fieldnames(BoundaryConditions))
		for (field,i) ∈ zip(fields, 1:length(fields)-1)
			BC_array[i] = haskey(bc_dict, field) ? bc_dict[field] : nothing
		end
		BoundaryConditions(BC_array..., zeros(4))
	end

	function evaluate(f::Vector{Function}, t::Union{Float64,Int64})
		return [f[i](t) for i ∈ 1:length(f)]
	end

	function solve_st(sys::System) #Stationary solver
		u = sys.S\evaluate(sys.f, 0)
		return x -> CubicHermiteSpline.CubicHermiteSplineInterpolation(
		    sys.x,u[1:2:end],u[2:2:end])(x), u
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
		
		u = zeros(nₓ,nₜ); u̇ = zeros(nₓ,nₜ); ü	= zeros(nₓ,nₜ);
		u[:,1] = IC[:,1]
		u̇[:,1] = IC[:,2]
		ü[:,1] = IC[:,3]
		for j=1:nₜ-1
			hⱼ = times[j+1] - times[j]
			uⱼ_star = u[:,j] + u̇[:,j]*hⱼ + (1/2 - β)*ü[:,j]*hⱼ^2
			u̇ⱼ_star = u̇[:,j] + (1 - γ)*ü[:,j]*hⱼ

			rhs = evaluate(sys.f, times[j+1]) - sys.S*uⱼ_star
			rhs[sys.par.BCs.indices] .= 0 # FIX: potential error here?

			ü[:,j+1] = (sys.M+β*hⱼ^2*sys.S)\rhs
			u̇[:,j+1] = u̇ⱼ_star + γ*ü[:,j+1]*hⱼ
			u[:,j+1] = uⱼ_star + β*ü[:,j+1]*hⱼ^2			
		end
		
		return [x -> CubicHermiteSpline.CubicHermiteSplineInterpolation(
			sys.x, u[1:2:end,j], u[2:2:end,j])(x) for j ∈ 1:length(times)], u
	end

	# Applies Dirichlet and Neumann BCs at x=0 and x=L
	function set_BCs(
					S::SparseArrays.SparseMatrixCSC{Float64,Int64}, 
					M::SparseArrays.SparseMatrixCSC{Float64,Int64}, 
					f::Vector{Function}, 
					N_u::Integer, 
					par::Parameters,
					x::Vector{Float64})

		cnt = 1
		display(par.BCs)
		if par.BCs.x_0 !== nothing
			S[N_u + 1, 1] = 1
			f[N_u + 1] = t -> par.BCs.x_0
			par.BCs.indices[cnt] = N_u + 1; cnt += 1
		end
		if par.BCs.xprime_0 !== nothing
			S[N_u + 2, 2] = 1
			f[N_u + 2] = t -> par.BCs.xprime_0
			par.BCs.indices[cnt] = N_u + 2; cnt += 1
		end
		if par.BCs.x_L !== nothing
			S[N_u + 3, end-1] = 1
			f[N_u + 3] = t -> par.BCs.x_L
			par.BCs.indices[cnt] = N_u + 3; cnt += 1
		end
		if par.BCs.xprime_L !== nothing
			S[N_u + 4, end] = 1
			f[N_u + 4] = t -> par.BCs.xprime_L
			par.BCs.indices[cnt] = N_u + 4; cnt += 1
		end

		h_0 = x[2]-x[1]
		h_L = x[end]-x[end-1]
		L = x[end]

		if par.BCs.Q_0 !== nothing
			S[N_u + 5, 1:4] = Q_at_point(0.0, h_0, par)
			f[N_u + 5] = t -> -par.BCs.Q_0 
			par.BCs.indices[cnt] = N_u + 5; cnt += 1
		end
		if par.BCs.M_0 !== nothing
			S[N_u + 6, 1:4] = -[6/h_0^2	 4/h_0	-6/h_0^2	 2/h_0] * par.EI(0)
			f[N_u + 6] = t -> -par.BCs.M_0 
			par.BCs.indices[cnt] = N_u + 6; cnt += 1
		end

		if par.BCs.Q_L !== nothing
			S[N_u + 7, end-3:end] = -Q_at_point(L, h_L, par)
			f[N_u + 7] = t -> -par.BCs.Q_L 
			par.BCs.indices[cnt] = N_u + 7; cnt += 1
		end
		if par.BCs.M_L !== nothing
			S[N_u + 8, end-3:end] = [6/h_L^2	 2/h_L	-6/h_L^2	 4/h_L] * par.EI(L)
			f[N_u + 8] = t -> par.BCs.M_L 
			par.BCs.indices[cnt] = N_u + 8; cnt += 1
		end

	# M[par.BCs.indices, :] .= 0 # FIX: Don't think this is correct assignment
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
		f_loc(h,p0) = [t -> GQ3(h,p->phi_0(h,p)[i]*par.q(p+p0,t)) for i ∈ 1:4]
		
		#Global Variables
		N_v = length(x) # Number of vertices
		N_e = N_v-1     # Number of elements
		N_u = N_v*2     # Number of unknowns
		N_bc = 8		# Number of (possible) Boundary Conditions

		#Global System
		M = SparseArrays.spzeros(Float64,N_u + N_bc,N_u)
		S = SparseArrays.spzeros(Float64,N_u + N_bc,N_u)
		f = Vector{Function}(undef, N_u + N_bc)
		for i ∈ 1:N_u + N_bc
			f[i] = t -> 0
		end

		#Element contributions
		for k ∈ 1:N_e
			i       = i_loc.+2*(k-1)
			h       = x[k+1]-x[k]
			
			M[i,i] += M_loc(h,x[k])
			S[i,i] += S_loc(h,x[k])
			for j ∈ 1:4
				i_f = j+2*(k-1)
				ftemp = f[i_f]
				f[i_f] = t -> ftemp(t) + f_loc(h,x[k])[j](t)
			end
		end

		# S embeds boundary conditions in an 8xN block at the bottom of the matrix.
		# First four rows are for x and x' BCs.
		# #Boundary Conditions
		set_BCs(S, M, f, N_u, par, x)

		#Packaging
		return System(par,x,S,M,f)
	end

	# Apply Q at a given point . Usually 0 or L for BCs
	# Returns a vector to place in our S matrix
	function Q_at_point(point::Float64, h::Float64, par::Parameters)
		EIprime = ForwardDiff.derivative(par.EI, point)

		S_vector = [6/h^2,	 4/h,	-6/h^2,	 2/h] * EIprime
		S_vector += [12/h^3,	 6/h^2,	-12/h^3,	 6/h^2] * par.EI(point)

		return S_vector
	end
end