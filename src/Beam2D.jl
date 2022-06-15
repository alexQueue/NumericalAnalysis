module Beam2D
	using SparseArrays,Printf #Stdlib imports
	import IterTools, Arpack, CubicHermiteSpline #External imports

    function problem_constructor_from_file(file::String)
        io = open(file, "r")
        txt = read(io, String)
        fl_as_list = split(txt,"\n")
        fl_as_list = strip.(fl_as_list) # Remove trailing whitespace
        filter!(x->x[1] != '#', fl_as_list) # Remove comments
        fl_as_list = String.(fl_as_list)
        fl_as_list = striplinecomment.(fl_as_list)
        indices = findall(x->x in ["NODES","CONNECTION","TYPE"], fl_as_list) 
        splitted = getindex.(Ref(fl_as_list), UnitRange.([1; indices .+ 1], [indices .- 1; length(fl_as_list)])) # Split into 3 lists
        nodes,connections,types = splitted[2:end] # 2:end because 1st element is empty as we split on "NODES"
    end

    function striplinecomment(a::String, cchars::String="#")
        b = strip(a)
        0 < length(cchars) || return b
        for c in cchars
            r = Regex(@sprintf "\\%c.*" c)
            b = replace(b, r => "")
        end
        strip(b)
    end

	# struct Problem
	# 	parameters ::NamedTuple{(:mu,:EI,:q),NTuple{3,Function}}
	# 	BCs        ::Dict{Tuple{Bool,Char},Float64} #(side,type)=>value
	# 	grid       ::Vector{Float64}

	# end

	# struct System
	# 	problem ::Problem
	# 	shape   ::Tuple{Int,Int}
	# 	Me      ::SparseMatrixCSC{Float64,Int64}
	# 	Se      ::SparseMatrixCSC{Float64,Int64}
	# 	qe      ::Vector{Float64}

	# 	function System(problem::Problem)
	# 		#Index calculations
	# 		n_w = 2*length(problem.grid) #Number of internal variables
	# 		n_b = 4                      #Number of boundary variables
	# 		n_c = length(problem.BCs)    #Number of boundary conditions
	# 		n_s = n_w + n_c              #Number of system equations
	# 		n_x = n_w + n_b              #Number of state variables

	# 		#Global FEM system
	# 		Me = spzeros(n_s,n_x)
	# 		Se = spzeros(n_s,n_x)
	# 		qe = zeros(Float64,n_s)
			
	# 		M = view(Me,1:n_w,1:n_w)     #Mass matrix
	# 		S = view(Se,1:n_w,1:n_w)     #Stiffness matrix
	# 		B = view(Se,1:n_w,n_w+1:n_x) #Boundary term matrix
	# 		C = view(Se,n_w+1:n_s,1:n_x) #Condition matrix
	# 		q = view(qe,1:n_w)           #Load vector
	# 		c = view(qe,n_w+1:n_s)       #Condition vector

	# 		#Physics
	# 		#Shape functions and derivatives on element [0,h], with p in [0,1]
	# 		phi_0(h,p) = [ p^2*(2*p-3)+1, p*(p-1)^2*h  ,-p^2*(2*p-3)  , p^2*(p-1)*h   ]
	# 		phi_1(h,p) = [ 6*p*(p-1)/h  , (p-1)*(3*p-1),-6*p*(p-1)/h  , p*(3*p-2)     ]
	# 		phi_2(h,p) = [ (12*p-6)/h^2 , (6*p-4)/h    ,-(12*p-6)/h^2 , (6*p-2)/h     ]
	# 		phi_3(h,p) = [ 12/h^3       , 6/h^2        ,-12/h^3       , 6/h^2         ]

	# 		#1D, 3 point Gaussian quadrature on element [0,h], with p in [0,1]
	# 		GQ3(h,f) = 5*h/18*f(1/2-sqrt(3/20)) +
	# 		           4*h/9 *f(1/2           ) +
	# 		           5*h/18*f(1/2+sqrt(3/20))
			
	# 		#Local System for element [o,o+h], with p in [0,1]
	# 		M_loc(h,o) = GQ3(h,p -> phi_0(h,p)*phi_0(h,p)'*problem.parameters.mu(o+h*p))
	# 		S_loc(h,o) = GQ3(h,p -> phi_2(h,p)*phi_2(h,p)'*problem.parameters.EI(o+h*p))
	# 		B_loc(h,o) = [-phi_1(h,0) -phi_0(h,0)  phi_1(h,1)  phi_0(h,1)] .* 
	# 		             (o .== problem.grid[[1,1,end-1,end-1]]')
	# 		q_loc(h,o) = GQ3(h,p -> phi_0(h,p)*problem.parameters.q(o+h*p))
			
	# 		#Contributions for each element
	# 		for (i,h,o) in zip(collect.(IterTools.partition(1:n_w,4,2)),diff(problem.grid),problem.grid)
	# 			M[i,i] += M_loc(h,o)
	# 			S[i,i] += S_loc(h,o)
	# 			B[i,:] += B_loc(h,o)
	# 			q[i]   += q_loc(h,o)
	# 		end

	# 		#Boundary conditions
	# 		for (i,((side,type),val)) in enumerate(problem.BCs)
	# 			j = type in ['H','G'] ?
	# 			    1     + (type=='G') + side*(n_w-2) :
	# 			    n_w+1 + (type=='Q') + side*2

	# 			C[i,j] = 1
	# 			c[i]   = val
	# 		end

	# 		#Packaging
	# 		return new(problem,(n_s,n_x),Me,Se,qe)
	# 	end
	# end

	# function u_to_Vh(grid::Vector{Float64},u::AbstractVector{Float64}) #Convert coefficients to Vh function
	# 	@assert length(u) == 2*length(grid) "Wrong number of coefficients for given grid"
		
	# 	return x -> CubicHermiteSpline.CubicHermiteSplineInterpolation(grid,u[1:2:end],u[2:2:end])(x)
	# end

	# function solve_st(sys::System) #Stationary solver
	# 	return u_to_Vh(sys.problem.grid,(sys.Se\sys.qe)[1:end-4])
	# end

	# function solve_tr_Newmark(sys::System,IC::Matrix{Float64},times::Vector{Float64})
	# 	@assert size(IC) == (sys.shape[2],3) "Wrong IC size for given system"
	# 	@assert length(times) >= 2 "Must have an initial and final time"
	# 	@assert all(diff(times) .> 0) "Times must be ascending"
		
	# 	n_t = length(times)
		
	# 	u_0 = Array{Float64,2}(undef,sys.shape[2],n_t)
	# 	u_1 = Array{Float64,2}(undef,sys.shape[2],n_t)
	# 	u_2 = Array{Float64,2}(undef,sys.shape[2],n_t)
		
	# 	u_0[:,1], u_1[:,1], u_2[:,1] = eachcol(IC)

	# 	beta  = 1/4
	# 	gamma = 1/2

	# 	for (t,h) in enumerate(diff(times))
	# 		u_0s = u_0[:,t] + h*u_1[:,t] + (0.5-beta)*h^2*u_2[:,t]
	# 		u_1s = u_1[:,t] + (1-gamma)*h*u_2[:,t]

	# 		u_2[:,t+1] = (sys.Me+beta*h^2*sys.Se)\(sys.qe-sys.Se*u_0s)
	# 		u_1[:,t+1] = u_1s + gamma*h*u_2[:,t+1]
	# 		u_0[:,t+1] = u_0s + beta*h^2*u_2[:,t+1]
	# 	end
		
	# 	return [u_to_Vh(sys.problem.grid,u) for u in eachcol(u_0[1:end-4,:])]
	# end

	# function get_numerical_eigenvalues(sys::System)

	# 	eigenvals, eigenvecs = Arpack.eigs(sys.Me, sys.Se)	#could be changed for SM
	# 	#may resort eigenvalues and eigenvecs
	# 	eigenval_no_zeros = eigenvals[ eigenvals .!= 0]
	# 	eigenvecs_no_zeros = eigenvecs[: , eigenvals .!= 0]
	
	# 	w_ks = 1 ./ sqrt.(eigenval_no_zeros)
	# 	w = (t, a, b) -> 0
	# 	for w_j in w_ks
	# 		w_i = (t, a, b) -> a .*cos(w_j*t) .+ (b/w_j).* sin(w_j*t)
	# 		w_temp = w
	# 		w(t, a, b) = w_temp(t, a, b) .+ w_i(t, a, b)
	# 	end

	# 	return real.(eigenval_no_zeros), real.(eigenvecs_no_zeros), w
	# end

	# function get_analytic_eigenvalues(system::System, pars, L)
	# 	n = 4
	# 	A = L^(-1/2)
	# 	x_j = [(x-0.5)*pi/L for x in 1:n]
	# 	k = [x_j[x]/L for x in 1:n ]
	# 	# analytic case: Assumption that EI and mu are constant
	# 	freq = (pars.EI(0)./pars.mu(0)) .* k.^4
	# 	w_j(x) =  [  
	# 		x -> A .* (
	# 			(cosh(k[j].*x) .- cos(k[j] .*x)) 
	# 			.- ((cosh(x_j[j]) .- cos(x_j[j]))./(sinh(x_j[j]) .+ sin(x_j[j])) 
	# 			.* (sinh(k[j] .*x) .- sin(k[j] .*x)))
	# 			) for j in 1:n
	# 		]

	# 	return freq, w_j
	# end

end