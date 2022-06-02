module Beam1D
	import SparseArrays, CubicHermiteSpline

	struct Problem
		parameters::NamedTuple{(:mu,:EI,:q),Tuple{Function,Function,Function}}
		BCs::Dict{Tuple{Bool,Char},Float64} #(side,type)=>value
		grid::Vector{Float64}

		Problem(p,b,g) = length(b) != 4 ?
		                 throw("Invalid BCs") :
		                 new(p,b,g)
	end

	struct System
		problem::Problem
		shape::Tuple{Int,Int}
		A::SparseArrays.SparseMatrixCSC{Float64,Int64}
		E::SparseArrays.SparseMatrixCSC{Float64,Int64}
		f::Vector{Float64}
	end

	function u_to_Vh(sys::System,u) #Convert coefficients to Vh function
		return x -> CubicHermiteSpline.CubicHermiteSplineInterpolation(
		            sys.problem.grid,u[1:2:end-4],u[2:2:end-4])(x)
	end

	function solve_st(sys::System) #Stationary solver
    u = sys.E\sys.f
		return u_to_Vh(sys,u), u
	end

	function solve_tr(sys::System,IC::Matrix{Float64},times::Vector{Float64})
		N_u = size(IC)[1]
		N_t = length(times)
		
		u_0 = Array{Float64,2}(undef,N_u,N_t)
		u_1 = Array{Float64,2}(undef,N_u,N_t)
		u_2 = Array{Float64,2}(undef,N_u,N_t)
		
		u_0[:,1] = IC[:,1]; u_1[:,1] = IC[:,2]; u_2[:,1] = IC[:,3]

		beta  = 1/4
		gamma = 1/2

		for (t,h) in enumerate(diff(times))
			u_0s = u_0[:,t] + h*u_1[:,t] + (0.5-beta)*h^2*u_2[:,t]
			u_1s = u_1[:,t] + (1-gamma)*h*u_2[:,t]

			u_2[:,t+1] = (sys.A+beta*h^2*sys.E)\(sys.f-sys.E*u_0s)
			u_1[:,t+1] = u_1s + gamma*h*u_2[:,t+1]
			u_0[:,t+1] = u_0s + beta*h^2*u_2[:,t+1]
		end
		
		return [u_to_Vh(sys,u) for u in eachcol(u_0)], u_0
	end

	function build(problem::Problem)
		#Shape functions and derivatives on element [0,h], with p in [0,1]
		phi_0(h,p) = [ p^2*(2*p-3)+1  p*(p-1)^2*h   -p^2*(2*p-3)    p^2*(p-1)*h   ]
		phi_1(h,p) = [ 6*p*(p-1)/h    (p-1)*(3*p-1) -6*p*(p-1)/h    p*(3*p-2)     ]
		phi_2(h,p) = [ (12*p-6)/h^2   (6*p-4)/h     -(12*p-6)/h^2   (6*p-2)/h     ]
		phi_3(h,p) = [ 12/h^3         6/h^2         -12/h^3         6/h^2         ]

		#1D, 3 point Gaussian quadrature on element [0,h], with p in [0,1]
		GQ3(h,f) = 5*h/18*f(1/2-sqrt(3/20)) +
		           4*h/9 *f(1/2           ) +
		           5*h/18*f(1/2+sqrt(3/20))

		#Local System for element [o,o+h], with p in [0,1]
		i_loc      = [1, 2, 3, 4]
		M_loc(h,o) = GQ3(h,p -> phi_0(h,p)'*phi_0(h,p)*problem.parameters.mu(o+h*p))
		S_loc(h,o) = GQ3(h,p -> phi_2(h,p)'*phi_2(h,p)*problem.parameters.EI(o+h*p))
		B_loc(h,o) = [-phi_1(h,0)' -phi_0(h,0)' phi_1(h,1)' phi_0(h,1)'] .* 
		             (o .== problem.grid[[1,1,end-1,end-1]]')
		q_loc(h,o) = GQ3(h,p -> phi_0(h,p)'*problem.parameters.q(o+h*p))

		#Global variables
		N_v = length(problem.grid) #Number of grid vertices
		N_e = N_v-1                #Number of elements
		N_w = N_v*2                #Number of internal unknowns
		N_b = 4                    #Number of boundary unknowns
		N_c = length(problem.BCs)  #Number of boundary conditions
		N_s = (N_w+N_c,N_w+N_b)    #System shape

		#Global system
		A = SparseArrays.spzeros(Float64,N_s[1],N_s[2])
		E = SparseArrays.spzeros(Float64,N_s[1],N_s[2])
		f = zeros(Float64,N_s[1])
		
		M = view(A,1:N_w,1:N_w)           #Mass matrix
		S = view(E,1:N_w,1:N_w)           #Stiffness matrix
		B = view(E,1:N_w,N_w+1:N_s[2])    #Boundary term matrix
		C = view(E,N_w+1:N_s[1],1:N_s[2]) #Condition matrix
		q = view(f,1:N_w)                 #Load vector
		c = view(f,N_w+1:N_s[1])          #Condition vector

		#Physics
		for k in 1:N_e
			i = i_loc.+2*(k-1)
			h = problem.grid[k+1] - problem.grid[k]

			M[i,i] += M_loc(h,problem.grid[k])
			S[i,i] += S_loc(h,problem.grid[k])
			B[i,:] += B_loc(h,problem.grid[k])
			q[i]   += q_loc(h,problem.grid[k])
		end

		#Boundary conditions
		for (i,((side,type),val)) in enumerate(problem.BCs)
			j = type in ['H','G'] ?
			    1     + (type=='G') + side*(N_w-2) :
			    N_w+1 + (type=='Q') + side*2

			C[i,j] = 1
			c[i]   = val
		end

		#Packaging
		return System(problem,N_s,A,E,f)
	end
end