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

		#Global System
		M = SparseArrays.spzeros(Float64,N_u,N_u)
		S = SparseArrays.spzeros(Float64,N_u,N_u)
		f = zeros(Float64,N_u)

		#Element contributions
		for k in 1:N_e
			i       = i_loc.+2*(k-1)
			h       = x[k+1]-x[k]
			
			M[i,i] += M_loc(h,x[k])
			S[i,i] += S_loc(h,x[k])
			f[i]   += f_loc(h,x[k])
		end

		#Boundary Conditions
		i      = [1,2,N_u-1,N_u] #Boundary indices
		S[i,:] = SparseArrays.sparse(LinearAlgebra.I,N_u,N_u)[i,:]
		f[i]   = par.BCs

		#Packaging
		return System(par,x,S,f)
	end
end
