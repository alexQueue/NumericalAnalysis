module Beam2D
	using SparseArrays,Printf,LinearAlgebra #Stdlib imports
	import IterTools, Arpack, CubicHermiteSpline #External imports

    mutable struct Node
        type::String
        coord::Vector{Float64}
        number::Int64

        force::Vector{Float64}
        moment::Float64
        movable_direction::Union{Vector{Float64},String}

        connecting_edges::Vector{Int64}

        function Node(
                    type::Union{String,SubString{String}},
                    coord::Vector{Float64},
                    number::Int64;
                    force::Vector{Float64} = [0.0,0.0],
                    moment::Float64 = 0.0,
                    movable_direction::Union{Vector{Float64},String} = "ALL"
                )
            new(type,coord,number,force,moment,movable_direction,[])
        end
    end

    struct Edge
        nodes::Vector{Node}
        grid::Vector{Float64}
        index_start::Int64
    end

    struct Problem
		EI::Float64
        nodes::Vector{Node}
		edges::Vector{Edge}
        shape::Tuple{Int,Int}

        function Problem(EI::Float64, nodes::Vector{Node}, edges::Vector{Edge})
            last_edge = edges[end]
            size = last_edge.index_start + length(last_edge.grid)*3 - 1
            shape = (size,size)

            for (edge,i) in zip(edges,1:length(edges))
                for node in edge.nodes
                    push!(node.connecting_edges, i)
                end
            end

            new(EI,nodes,edges,shape)
        end
	end

    function problem_constructor(EI::Float64, file::String)
        # Read file and save initial lists
        io = open(file, "r")
        txt = read(io, String)
        fl_as_list = split(txt,"\n")
        fl_as_list = strip.(fl_as_list) # Remove trailing whitespace
        filter!(x->x[1] != '#', fl_as_list) # Remove comments
        fl_as_list = String.(fl_as_list)
        fl_as_list = striplinecomment.(fl_as_list)
        indices = findall(x->x in ["NODES","EDGES","TYPE"], fl_as_list) 
        splitted = getindex.(Ref(fl_as_list), UnitRange.([1; indices .+ 1], [indices .- 1; length(fl_as_list)])) # Split into 3 lists
        nodes,edges,types = splitted[2:end] # 2:end because 1st element is empty as we split on "NODES"

        nodes = strlist_to_type(Float64, nodes)
        edges = strlist_to_type(Int, edges)
        
        types = get_node_type_list(types)
        # Adjacency = create_adjacency_matrix(nodes, edges)

        Nodes = Vector{Node}(undef, length(nodes))
        for (node,type,i) in zip(nodes,types,1:length(nodes))
            if type[2] == "FORCE"
                force = parse.(Float64, split(type[3]))
                moment = parse.(Float64, type[4])
                Nodes[i] = Node(type[2], node, i, force=force, moment=moment)
            elseif type[2] == "MOVABLE"
                movable_direction = parse.(Float64, split(type[3]))
                Nodes[i] = Node(type[2], node, i, movable_direction=movable_direction)
            else
                Nodes[i] = Node(type[2], node, i)
            end
        end

        index_cnt = 1
        Edges = Vector{Edge}(undef, length(edges))
        for (edge,i) in zip(edges,1:length(edges))
            grid = [0.0, norm(Nodes[edge[1]].coord - Nodes[edge[2]].coord)] # 0 -> length of beam
            Edges[i] = Edge([ Nodes[edge[1]], Nodes[edge[2]] ], grid, index_cnt)
            index_cnt += length(grid)*3 # v_1 -> v_n, and w_1 -> w_2n makes 3n
        end

        Problem(EI,Nodes,Edges)
    end
    
    function get_node_type_list(types::Vector)
        force_re = r" +\[(.*?)\] *?| "
        forces = []
        # Get the prescribed forces/momenta/directed movement 
        # from the types list
        for type in types
            matches = eachmatch(force_re, type, overlap = true)
            for match in matches
                if match[1] !== nothing
                    push!(forces, match.captures) 
                end
            end
        end
        
        types = [split(x, force_re) for x in types]

        # Add the correct forces/momenta/directed movement 
        # to the list of type of each node
        for (list,i) in zip(types,1:length(types))
            if list[2] == "FORCE"
                types[i][3] = forces[1][1]
                types[i][4] = forces[2][1]
                popfirst!(forces)
                popfirst!(forces)
            elseif list[2] == "MOVABLE"
                types[i][3] = forces[1][1]
                popfirst!(forces)
            end
        end
        types
    end

    function fixed_bearings_nodes(types)
        fb_nodes = 0
        for node in types
        end
    end

    function strlist_to_type(type::Type, data::Vector)
        data = split.(data, " ")
        [parse.(type, x) for x in data]
    end

    function create_adjacency_matrix(nodes::Vector, edges::Vector)
        size = length(nodes)
        Adjacency = spzeros(size,size)
        for edge in edges
            d = norm(nodes[edge[1]]-nodes[edge[2]])
            Adjacency[edge[1],edge[2]] = d
        end
        Adjacency + transpose(Adjacency)
    end

    # Copied from internet obviously
    function striplinecomment(a::String, cchars::String="#")
        b = strip(a)
        0 < length(cchars) || return b
        for c in cchars
            r = Regex(@sprintf "\\%c.*" c)
            b = replace(b, r => "")
        end
        strip(b)
    end

    function C_matrix_construction(Problem)
        # Count size of cᵀ, denoted r
        r = length(Problem.edges)*3
        r = 0
        for node in Problem.nodes
            if node.type == "FIXED"
                r += 3
            elseif node.type == "MOVABLE"
                n_cnct_edges = length(node.connecting_edges)
                r += n_cnct_edges # one bearing condition per connecting edge
                if n_cnct_edges >= 2
                    r += 0
                end
            else
                r += (length(node.connecting_edges) - 1)*3
            end
        end

        C = spzeros(r,Problem.shape[1])
        i = 1
        for node in Problem.nodes
            if node.type == "FIXED"
                j = Problem.edges[node.connecting_edges[1]].index_start
                C[i,j] = 1
                i += 1

                j += length(Problem.edges[node.connecting_edges[1]].grid)
                C[i,j] = 1
                i += 1
                
                j += 1
                C[i,j] = 1
                i += 1
            elseif node.type == "MOVABLE"
                for edge in Problem.edges[node.connecting_edges]
                    phi = edge_angle(edge)
                    order = edge_node_order(edge, node)

                    # No movement in movable direction means the dot product of movement
                    # vector with flipped displacement vector is 0
                    # Flipped because displacement vector should be same direction as movement
                    j = linking_index(edge, order)
                    dx = node.movable_direction[1]; dy = node.movable_direction[2]
                    
                    C[i,j] = [dx*cos(phi) + dy*sin(phi), -dx*sin(phi) + dy*cos(phi)]
                    i += 1
                end
            else # "FORCE/FREE"
                i = connecting_edges_conditions!(Problem, node, C, i)
            end
        end
        C
    end

    function edge_node_order(edge, node)
        order = edge.nodes[1].number == node.number ? "last" : "first"
    end

    function connecting_edges_conditions!(Problem, node, C, i)
        first_edge = Problem.edges[node.connecting_edges[1]]
        phi_1 = edge_angle(first_edge)
        # Find out if connecting edge ends or begins at the node
        first_order = edge_node_order(first_edge, node)
        
        for rem_edge in Problem.edges[node.connecting_edges[2:end]]
            rem_order = edge_node_order(rem_edge, node)
            
            # Stiffness condition
            j1 = stiffness_index(first_edge, first_order)
            j2 = stiffness_index(rem_edge, rem_order)
            C[i,[j1,j2]] = [1,-1]
            i += 1

            # Linking condition
            phi_2 = edge_angle(rem_edge)
            j1 = linking_index(first_edge, first_order)
            j2 = linking_index(rem_edge, rem_order)
            
            C[i,[j1...,j2...]] = [cos(phi_1), -sin(phi_1), -cos(phi_2), sin(phi_2)]
            i += 1
            
            C[i,[j1...,j2...]] = [sin(phi_1), cos(phi_1), -sin(phi_2), -cos(phi_2)]
            i += 1
        end
        return i
    end

    # Returns the indices for the linking condition for the edge with order "last" or "first"
    function linking_index(edge, order)
        i = order == "first" ? 
            [edge.index_start, 
                edge.index_start + length(edge.grid)] :
            [edge.index_start + length(edge.grid) - 1,
                edge.index_start + length(edge.grid)*3 - 2]
        return i
    end

    function stiffness_index(edge, order)
        i = order == "first" ? 
            edge.index_start + length(edge.grid)*3 - 1 :
            edge.index_start + length(edge.grid) + 1 
        return i
    end

    function edge_angle(edge)
        dy = edge.nodes[2].coord[2] - edge.nodes[1].coord[2]
        dx = edge.nodes[2].coord[1] - edge.nodes[1].coord[1]
        atan(dy,dx)
    end

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