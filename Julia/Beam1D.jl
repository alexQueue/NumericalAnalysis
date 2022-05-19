module Beam1D
	import SparseArrays, LinearAlgebra, CubicHermiteSpline

	mutable struct Parameters
		mu::Function
		EI::Function
		q::Function
		BCs::Vector{Float64}
	end

	struct System
		par::Parameters
		x::Vector{Float64}
		S::SparseArrays.SparseMatrixCSC{Float64,Int64}
		M::SparseArrays.SparseMatrixCSC{Float64,Int64}
		f::Vector{Function}
	end

	function evaluate(f,t)
		return [f[i](t) for i ∈ 1:length(f)]
	end

	function solve_st(sys::System) #Stationary solver
		u = sys.S\sys.f(0)
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

		u = zeros(nₓ,nₜ); u̇ = zeros(nₓ,nₜ); ü	= zeros(nₓ,nₜ);
		u[:,1] = IC[:,1]
		u̇[:,1] = IC[:,2]
		ü[:,1] = IC[:,3]
		for j=1:nₜ-1
			hⱼ = times[j+1] - times[j]
			uⱼ_star = u[:,j] + u̇[:,j]*hⱼ + (1/2 - β)*ü[:,j]*hⱼ^2
			u̇ⱼ_star = u̇[:,j] + (1 - γ)*ü[:,j]*hⱼ
	
			rhs = evaluate(sys.f, times[j+1]) - sys.S*uⱼ_star
			N_u = length(sys.x)*2
			i = [1,2,N_u,N_u-1]
			rhs[i] .= 0

			ü[:,j+1] = (sys.M+β*hⱼ^2*sys.S)\rhs
			u̇[:,j+1] = u̇ⱼ_star + γ*ü[:,j+1]*hⱼ
			u[:,j+1] = uⱼ_star + β*ü[:,j+1]*hⱼ^2
		end
		return [x -> CubicHermiteSpline.CubicHermiteSplineInterpolation(
			sys.x, u[1:2:end,j], u[2:2:end,j])(x) for j ∈ 1:length(times)]
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
		N_v = length(x) #Number of vertices
		N_e = N_v-1     #Number of elements
		N_u = N_v*2     #Number of unknowns

		#Global System
		M = SparseArrays.spzeros(Float64,N_u,N_u)
		S = SparseArrays.spzeros(Float64,N_u,N_u)
		f = Vector{Function}(undef, N_u)
		for i ∈ 1:N_u
			f[i] = t->0
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

		#Boundary Conditions
		i      = [1,2,N_u-1,N_u] #Boundary indices
		S[i,:] = SparseArrays.sparse(LinearAlgebra.I,N_u,N_u)[i,:]
		M[i,:] .= 0
		for (j,k) ∈ zip(i,length(i))
			f[j]   = t->par.BCs[k]
		end

		#Packaging
		return System(par,x,S,M,f)
	end
end