# Beam1D.jl
# Authors: Alex Quinlan, Pia Callmer, Nicola Sabbadini & Henry Jacobson
#
# TU Berlin
# Project Numerical Analysis SoSe 2022

module Beam1D
	using SparseArrays, LinearAlgebra #Stdlib imports
	import IterTools, Arpack #External imports

    """
        struct Problem
    
    The parameters (mu, EI and q(x)) defined for some specific setup, as well
    as the boundary conditions defined as a dictionary with keys being
    (Bool,Char), defining the side and type, with the values being the prescribed
    value there. The grid is the 1D points the beam is defined on.
    """
	struct Problem
		parameters ::NamedTuple{(:mu, :EI, :q), NTuple{3,Function}}
		BCs        ::Dict{Tuple{Bool,Char}, Float64} #(side, type)=>value
		grid       ::Vector{Float64}

		function Problem(parameters,BCs,grid)
			p = new(parameters,BCs,grid)

			@assert length(p.BCs) == 4 "Must have 4 BCs"

			return p
		end
	end

    """
        struct Problem

    The shape of the matrices, as well as the Me, Se and qe
    matrices.

    Me => Mass matrix
    Se => Stiffness matrix
    qe => Force vector
    """
	struct System
		problem ::Problem
		shape   ::Tuple{Int,Int}
		Me      ::SparseMatrixCSC{Float64,Int64}
		Se      ::SparseMatrixCSC{Float64,Int64}
		qe      ::Vector{Float64}

		function System(problem::Problem)
			#Index calculations
			n_w = 2*length(problem.grid) #Number of internal variables
			n_b = 4                      #Number of boundary variables
			n_c = length(problem.BCs)    #Number of boundary conditions
			n_s = n_w + n_c              #Number of system equations
			n_x = n_w + n_b              #Number of state variables

			#Global FEM system
			Me = spzeros(n_s,n_x)
			Se = spzeros(n_s,n_x)
			qe = zeros(Float64,n_s)
			
			M = view(Me,1:n_w,1:n_w)     #Mass matrix
			S = view(Se,1:n_w,1:n_w)     #Stiffness matrix
			B = view(Se,1:n_w,n_w+1:n_x) #Boundary term matrix
			C = view(Se,n_w+1:n_s,1:n_x) #Condition matrix
			q = view(qe,1:n_w)           #Load vector
			c = view(qe,n_w+1:n_s)       #Condition vector

			#Physics
			#Shape functions and derivatives on element [0,h], with p in [0,1]
			phi_0(h,p) = [ p^2*(2*p-3)+1, p*(p-1)^2*h  ,-p^2*(2*p-3)  , p^2*(p-1)*h   ]
			phi_1(h,p) = [ 6*p*(p-1)/h  , (p-1)*(3*p-1),-6*p*(p-1)/h  , p*(3*p-2)     ]
			phi_2(h,p) = [ (12*p-6)/h^2 , (6*p-4)/h    ,-(12*p-6)/h^2 , (6*p-2)/h     ]
			phi_3(h,p) = [ 12/h^3       , 6/h^2        ,-12/h^3       , 6/h^2         ]

			#1D, 3 point Gaussian quadrature on element [0,h], with p in [0,1]
			GQ3(h,f) = 5*h/18*f(1/2-sqrt(3/20)) +
			           4*h/9 *f(1/2           ) +
			           5*h/18*f(1/2+sqrt(3/20))
			
			#Local System for element [o,o+h], with p in [0,1]
			M_loc(h,o) = GQ3(h,p -> phi_0(h,p)*phi_0(h,p)'*problem.parameters.mu(o+h*p))
			S_loc(h,o) = GQ3(h,p -> phi_2(h,p)*phi_2(h,p)'*problem.parameters.EI(o+h*p))
			B_loc(h,o) = [-phi_1(h,0) -phi_0(h,0)  phi_1(h,1)  phi_0(h,1)] .* 
			             (o .== problem.grid[[1,1,end-1,end-1]]')
			q_loc(h,o) = GQ3(h,p -> phi_0(h,p)*problem.parameters.q(o+h*p))
			
			#Contributions for each element
			for (i,h,o) in zip(collect.(IterTools.partition(1:n_w,4,2)), diff(problem.grid), problem.grid)
				M[i,i] += M_loc(h,o)
				S[i,i] += S_loc(h,o)
				B[i,:] += B_loc(h,o)
				q[i]   += q_loc(h,o)
			end

			#Boundary conditions
			for (i,((side,type),val)) in enumerate(problem.BCs)
				j = type in ['H','G'] ?
				    1     + (type=='G') + side*(n_w-2) :
				    n_w+1 + (type=='Q') + side*2

				C[i,j] = 1
				c[i]   = val
			end

			#Packaging
			return new(problem,(n_s,n_x),Me,Se,qe)
		end
	end

    """
        function u_to_Vh(grid::Vector{Float64}, u::AbstractVector{Float64})

    
    Converts a solution `u` to parametric curves using the solution `u` and the
    specified `grid`.

    Returns two vectors of the x-functions and the y-functions, parametrized
    from 0 to 1, for each element in the grid.
    """
	function u_to_Vh(grid::Vector{Float64}, u::AbstractVector{Float64}) #Convert coefficients to Vh function
		@assert length(u) == 2*length(grid)+4 "Wrong number of coefficients for given grid"

		phi(t) = [t^2*(2*t-3)+1,t*(t-1)^2,-t^2*(2*t-3),t^2*(t-1)]

		xs = Vector{Function}(undef,length(grid)-1)
		ys = Vector{Function}(undef,length(grid)-1)

		for (i,(x0,x1),(y0,g0,y1,g1)) in zip(1:length(grid)-1,
		                                     IterTools.partition(grid      ,2,1),
		                                     IterTools.partition(u[1:end-4],4,2))
			h = x1-x0
			xs[i] = t -> dot([x0,   h,x1,   h],phi(t))
			ys[i] = t -> dot([y0,g0*h,y1,g1*h],phi(t))
		end

		return xs, ys
	end

    """
        function solve_st(sys::System)

    Solves the system `sys` and returns a vector of parametrized curves.
    """
	function solve_st(sys::System)
		return u_to_Vh(sys.problem.grid,(sys.Se\sys.qe))
	end

    """
        function solve_dy_Newmark(sys::System, IC::Matrix{Float64}, times::Vector{Float64})

    Calculates the solution of the system `sys` for specifiec initial conditions `IC`, containing
    the values for u, u_t and u_tt, and the discretized times for which the solution should
    be given. Returns a vector of vector of parametrized curves.
    """
	function solve_dy_Newmark(sys::System, IC::Matrix{Float64}, times::Vector{Float64})
		@assert size(IC) == (sys.shape[2],3) "Wrong IC size for given system"
		@assert length(times) >= 2 "Must have an initial and final time"
		@assert all(diff(times) .> 0) "Times must be ascending"
		
		n_t = length(times)
		
		u_0 = Array{Float64,2}(undef,sys.shape[2],n_t)
		u_1 = Array{Float64,2}(undef,sys.shape[2],n_t)
		u_2 = Array{Float64,2}(undef,sys.shape[2],n_t)
		
		u_0[:,1], u_1[:,1], u_2[:,1] = eachcol(IC)

		beta, gamma = (1/4,1/2)

		for (t,h) in enumerate(diff(times))
			u_0s = u_0[:,t] + h*u_1[:,t] + (0.5-beta)*h^2*u_2[:,t]
			u_1s = u_1[:,t] + (1-gamma)*h*u_2[:,t]

			u_2[:,t+1] = (sys.Me+beta*h^2*sys.Se)\(sys.qe-sys.Se*u_0s)
			u_1[:,t+1] = u_1s + gamma*h*u_2[:,t+1]
			u_0[:,t+1] = u_0s + beta*h^2*u_2[:,t+1]
		end
	
		return [u_to_Vh(sys.problem.grid,u) for u in eachcol(u_0)]
	end

    """
        function get_vibrations(sys::System,n_m::Int64=4)

    Calculates the vibrations of `sys` using the generalized eigenvalues,
    `n_m` determines how many eigenvalues to use.
    """
	function get_vibrations(sys::System, n_m::Int64=4)
		@warn "Boundary conditions and load assumed to be 0"

		evals, evecs = (-1,0)
		while any(evals .< 0)
			evals, evecs = real.(Arpack.eigs(sys.Me,sys.Se,nev=n_m))
		end

		freqs = evals.^(-0.5)
		modes = u_to_Vh.(Ref(sys.problem.grid),eachcol(evecs))
		
		return evals, evecs, freqs, modes
	end

    function get_eigenvalues_numerical(sys::System, n_m::Int64=4)
        evals, evecs, freqs, modes = get_vibrations(sys, n_m)

        w(t,a,b) = sum(a.*cos.(freqs.*t) .+ (b./freqs) .* sin.(freqs.*t))

        return evals, evecs, w
    end

    function get_eigenvalues_analytical(sys::System, pars, L, n_m)
        A = L^(-1/2)
        x_j = [(x-0.5)*pi/L for x in 1:n_m]
        k = [x_j[x]/L for x in 1:n_m]
        
        # analytic case: Assumption that EI and  mu are constant
        freq = pars.EI(0) / pars.mu(0) .* k.^4
        w_j(x) = [
            x -> A .* (
                (cosh(k[j].*x) .- cos(k[j] .*x)) 
                .- ((cosh(x_j[j]) .- cos(x_j[j]))./(sinh(x_j[j]) .+ sin(x_j[j])) 
                .* (sinh(k[j] .*x) .- sin(k[j] .*x)))
            ) for j in 1:n_m
        ]

        return freq, w_j 
    end

    """
        function solve_dy_eigen(sys::System, n_m::Int64=4)
    """
	function solve_dy_eigen(sys::System,n_m::Int64=4)
		evals, evecs, freqs, modes = get_vibrations(sys,n_m) 
		
		function get_sol(IC::Matrix{Float64})
			@assert size(IC) == (sys.shape[2],2) "Wrong IC size for given system"

			as = evecs\IC[:,1]
			bs = (evecs\IC[:,2])./freqs

			sol(t) = u_to_Vh(sys.problem.grid, evecs*(as.*cos.(freqs.*t).+bs.*sin.(freqs.*t)))

			return sol
		end

		return get_sol
	end

    """
        function eval(vec::Vector{Function}, x::Float64)
        
    Evaluate a vector of functions at x to be able to only use the midpoint of
    the xs,ys functions in a parametric curve for plotting.
    """
    function eval(vec::Vector{Function}, x::Float64)
        ret_vec = zeros(length(vec))
        for i in 1:length(vec)
            ret_vec[i] = vec[i](x)
        end
        return ret_vec
    end
end
