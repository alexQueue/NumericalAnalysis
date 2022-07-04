module Beam2D
	using SparseArrays,Printf,LinearAlgebra #Stdlib imports
	import IterTools, Arpack #External imports

    using PyCall

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
            if movable_direction != "ALL"
                movable_direction = movable_direction / norm(movable_direction)
            end
            new(type,coord,number,force,moment,movable_direction,[])
        end
    end

    struct Edge
        nodes::Vector{Node}
        grid::Vector{Float64}
        gridlen::Int64
        len::Float64
        index_start::Int64
    end

    struct Problem
		E::Function
		I::Function
        A::Function
        mu::Function
        nodes::Vector{Node}
		edges::Vector{Edge}
        shape::Tuple{Int,Int}

        function Problem(
                        E::Function, I::Function, A::Function, mu::Function,
                        nodes::Vector{Node}, edges::Vector{Edge}
                        )
            last_edge = edges[end]
            size = last_edge.index_start + length(last_edge.grid)*3 - 1
            shape = (size,size)

            for (i, edge) in enumerate(edges)
                for node in edge.nodes
                    push!(node.connecting_edges, i)
                end
            end

            new(E,I,A,mu,nodes,edges,shape)
        end
	end

    function problem_constructor(file::String)
        # Read file and save initial lists
        io = open(file, "r")
        txt = read(io, String)
        fl_as_list = split(txt,"\n")
        fl_as_list = strip.(fl_as_list) # Remove trailing whitespace
        filter!(x->x[1] != '#', fl_as_list) # Remove comments
        fl_as_list = String.(fl_as_list)
        fl_as_list = striplinecomment.(fl_as_list)
        indices = findall(x->x in ["NODES","EDGES","TYPE","PARAMETERS"], fl_as_list) 
        splitted = getindex.(Ref(fl_as_list), UnitRange.([1; indices .+ 1], [indices .- 1; length(fl_as_list)])) # Split into 4 lists
        nodes,edges,types,params = splitted[2:end] # 2:end because 1st element is empty as we split on "NODES"

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
        for (i,edge) in enumerate(edges)
            L = norm(Nodes[edge[1]].coord - Nodes[edge[2]].coord)
            gridpoints = 2
            grid = collect(LinRange(0.0,L,gridpoints)) # 0 -> length of beam with 
            Edges[i] = Edge([ Nodes[edge[1]], Nodes[edge[2]] ], grid, length(grid), L, index_cnt)
            index_cnt += length(grid)*3 # v_1 -> v_n, and w_1 -> w_2n makes 3n
        end

        parameters = Dict()
        for param in params
            str,val = split(param, " ", limit=2)
            parameters[str] = eval(Meta.parse("x -> " * val))
        end
        E = parameters["E"]; I = parameters["I"]; A = parameters["A"]; mu = parameters["mu"]

        Problem(E,I,A,mu,Nodes,Edges)
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
        for (i, list) in enumerate(types)
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

    function strlist_to_type(type::Type, data::Vector)
        data = split.(data, " ")
        [parse.(type, x) for x in data]
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

    # TODO - Add values to vector of values of constraints herein as well
    function C_matrix_construction(Problem)
        # Count size of cᵀ, denoted r
        r = 0
        for node in Problem.nodes
            if node.type == "FIXED"
                r += length(node.connecting_edges)*3 # 3 conditions per connecting edge
                # r += (length(node.connecting_edges) - 1)*3 # 3 conditions per pair of edges
            elseif node.type == "MOVABLE"
                n_cnct_edges = length(node.connecting_edges)
                r += n_cnct_edges # one bearing condition per connecting edge
                if n_cnct_edges >= 2
                    r += (n_cnct_edges - 1)*3
                end
            else
                r += (length(node.connecting_edges) - 1)*3 # 3 conditions per pair of edges
            end
        end

        C = spzeros(Problem.shape[1],r)
        f = spzeros(Problem.shape[1])

        i = 1
        for node in Problem.nodes
            if node.type == "FIXED"
                for edge in Problem.edges[node.connecting_edges]
                    j1,j2,j3 = fixed_index(edge, node)
                    
                    # v
                    C[j1,i] = 1
                    i += 1
                    
                    # w 
                    C[j2,i] = 1
                    i += 1
                    
                    # w'
                    C[j3,i] = 1
                    i += 1

                    # Add the connecting edges as well!
                    # i = connecting_edges_conditions!(Problem, node, C, i)
                end
            elseif node.type == "MOVABLE"
                for edge in Problem.edges[node.connecting_edges]
                    phi = edge_angle(edge)

                    # No movement in movable direction means the dot product of movement
                    # vector with flipped displacement vector is 0
                    # Flipped because displacement vector should be same direction as movement
                    j = linking_index(edge, node)
                    dx = node.movable_direction[1]; dy = node.movable_direction[2]
                    
                    C[j,i] = [dx*cos(phi) + dy*sin(phi), -dx*sin(phi) + dy*cos(phi)]
                    i += 1
                end
                # Even movable nodes can have more than one connecting edge 
                # and then we need stiffness and linking again
                # Q: Should stiffness condition hold here really?
                i = connecting_edges_conditions!(Problem, node, C, i)
            else # "FORCE/FREE"
                if node.type == "FORCE"
                    edge = Problem.edges[node.connecting_edges[1]]
                    angle = edge_angle(edge)
                    j1,j2 = linking_index(edge, node)
                    fx = node.force[1]; fy = node.force[2]
                    f[[j1,j2]] = [fx*cos(angle) + fy*sin(angle), -fx*sin(angle) + fy*cos(angle)]
                end
                i = connecting_edges_conditions!(Problem, node, C, i)
            end
        end
        C,f
    end

    function connecting_edges_conditions!(Problem, node, C, i)
        first_edge = Problem.edges[node.connecting_edges[1]]
        phi_1 = edge_angle(first_edge)
        
        for rem_edge in Problem.edges[node.connecting_edges[2:end]]
            # Stiffness condition
            j1 = stiffness_index(first_edge, node)
            j2 = stiffness_index(rem_edge, node)
            C[[j1,j2],i] = [1,-1]
            i += 1
            
            # Linking condition
            phi_2 = edge_angle(rem_edge)
            j1 = linking_index(first_edge, node)
            j2 = linking_index(rem_edge, node)
            
            C[[j1...,j2...],i] = [cos(phi_1), -sin(phi_1), -cos(phi_2), sin(phi_2)]
            i += 1
            
            C[[j1...,j2...],i] = [sin(phi_1), cos(phi_1), -sin(phi_2), -cos(phi_2)]
            i += 1
        end
        return i
    end

    function edge_node_order(edge, node)
        edge.nodes[1].number == node.number ? "first" : "last"
    end

    function fixed_index(edge, node)
        order = edge_node_order(edge, node)
        inds = order == "first" ? 
            [
                edge.index_start, 
                edge.index_start + edge.gridlen,
                edge.index_start + edge.gridlen+1
            ] : 
            [
                edge.index_start + edge.gridlen-1,
                edge.index_start + 2*edge.gridlen,
                edge.index_start + 2*edge.gridlen+1
            ]
        return inds
    end

    # Returns the indices for the linking condition for the edge with order "last" or "first"
    # Is the same indices for the force (i think)
    function linking_index(edge, node)
        order = edge_node_order(edge, node)
        inds = order == "first" ? 
            [edge.index_start, 
                edge.index_start + edge.gridlen] :
            [edge.index_start + edge.gridlen - 1,
                edge.index_start + edge.gridlen*3 - 2]
        return inds
    end

    function stiffness_index(edge, node)
        order = edge_node_order(edge, node)
        inds = order == "last" ? 
            edge.index_start + edge.gridlen*3 - 1 :
            edge.index_start + edge.gridlen + 1 
        return inds
    end

    function edge_angle(edge::Edge)
        dy = edge.nodes[2].coord[2] - edge.nodes[1].coord[2]
        dx = edge.nodes[2].coord[1] - edge.nodes[1].coord[1]
        atan(dy,dx)
    end

	struct System
		problem ::Problem
        shape   ::Tuple{Int,Int}
		Me      ::SparseMatrixCSC{Float64,Int64}
		Se      ::SparseMatrixCSC{Float64,Int64}
		qe      ::Vector{Float64}

		function System(problem::Problem)
            n = sum([edge.gridlen*3 for edge in problem.edges])

            C,f = C_matrix_construction(problem)
            r = size(C)[2]

            Me = spzeros(n_iv+r,n_iv+r)
            Se = spzeros(n_iv+r,n_iv+r)
            qe = zeros(Float64,n_iv+r)

            Se[1:n,n+1:n+r] = C
            Se[n+1:n+r,1:n] = Transpose(C)
            qe[1:n] = f

            # Physics
            # Shape functions and derivatives on element [0,h], with t in [0,1]
            phi_T_0(h,t) = [ t^2*(2*t-3)+1, t*(t-1)^2*h  ,-t^2*(2*t-3)  , t^2*(t-1)*h   ]
            phi_T_1(h,t) = [ 6*t*(t-1)/h  , (t-1)*(3*t-1),-6*t*(t-1)/h  , t*(3*t-2)     ]
            phi_T_2(h,t) = [ (12*t-6)/h^2 , (6*t-4)/h    ,-(12*t-6)/h^2 , (6*t-2)/h     ]
            phi_T_3(h,t) = [ 12/h^3       , 6/h^2        ,-12/h^3       , 6/h^2         ]
            phi_L_0(h,t) = [ t            , 1-t]
            phi_L_1(h,t) = [ 1/h          ,-1/h]

            # 1D, 3 point Gaussian quadrature on element [0,h], with t in [0,1]
			GQ3(h,f) = 5*h/18*f(1/2-sqrt(3/20)) +
                       4*h/9 *f(1/2           ) +
                       5*h/18*f(1/2+sqrt(3/20))

            EI(x) = problem.E(x) * problem.I(x)
            EA(x) = problem.E(x) * problem.A(x)

            # Longitudinal - Transversal ordering
            E_s = [EA,EI]
            M_base_fncs = [phi_L_0, phi_T_0]
            S_base_fncs = [phi_L_1, phi_T_2]

			# Local System for element [o,o+h], with t in [0,1]
			M_loc(h,o,base_fnc) = GQ3(h,t -> base_fnc(h,t)*base_fnc(h,t)'*problem.mu(o+h*t))
            S_loc(h,o,base_fnc,E_) = GQ3(h,t -> base_fnc(h,t)*base_fnc(h,t)'*E_(o+h*t))

            for edge in problem.edges
                partitions = [
                    # Longitudinal equation
                    edge.index_start:edge.index_start+edge.gridlen-1,
                    # Transversal equation
                    edge.index_start+edge.gridlen:edge.index_start+edge.gridlen*3-1,
                ]
                for j in 1:2 # Longitudinal / Transversal
                    partition = partitions[j]
                    Mbase_fnc = M_base_fncs[j]
                    Sbase_fnc = S_base_fncs[j]
                    E_ = E_s[j]

                    M = view(Me, partition, partition)
                    S = view(Se, partition, partition)
                    
                    for (i,h,o) in zip(collect.(IterTools.partition(1:n,2j,j)),diff(edge.grid),edge.grid)
                        M[i,i] += M_loc(h,o,Mbase_fnc)
                        S[i,i] += S_loc(h,o,Sbase_fnc,E_)
                    end
                end
            end
            shape = (n,r)

            new(problem,shape,Me,Se,qe)
        end
	end

    function getC(sys::System)
        n,r = sys.shape
        sys.Se[1:n,n+1:n+r]
    end

    function nspace(A::AbstractMatrix)
        sympy = pyimport("sympy")
        obj = sympy.Matrix(A).nullspace()
        convert.(Vector{Float64},obj)
    end

    function non0inds(v::Vector{Float64})
        B = []
        for (i,_) in enumerate(v)
            if v[i] != 0
                push!(B,i)
            end
        end
        B
    end

    function u_to_Vh(problem::Problem,u::AbstractVector{Float64}) #Convert coefficients to Vh function
        phi(t) = [t^2*(2*t-3)+1,t*(t-1)^2,-t^2*(2*t-3),t^2*(t-1)]

        xs = Vector{Function}(undef,length(problem.edges))
        ys = Vector{Function}(undef,length(problem.edges))

        for (i,edge) in enumerate(problem.edges)
            (x0,x1,y0,g0,y1,g1) = u[edge.index_start.+(0:5)]

            e = edge.nodes[2].coord - edge.nodes[1].coord
            R = [e[1] -e[2]; e[2] e[1]]./norm(e)
            h = norm(e) + x0 + x1

            q0_x,q0_y = edge.nodes[1].coord .+ R * [x0,y0]
            q1_x,q1_y = edge.nodes[2].coord .+ R * [x1,y1]
            u0_x,u0_y = R * [1,g0] .* h
            u1_x,u1_y = R * [1,g1] .* h

            xs[i] = t -> LinearAlgebra.dot([q0_x,u0_x,q1_x,u1_x],phi(t))
            ys[i] = t -> LinearAlgebra.dot([q0_y,u0_y,q1_y,u1_y],phi(t))
        end

        return xs, ys
    end

    function solve_st(sys::System) #Stationary solver
        return u_to_Vh(sys.problem.grid,(sys.Se\sys.qe))
    end

    function solve_dy_Newmark(sys::System,IC::Matrix{Float64},times::Vector{Float64})
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

    function get_vibrations(sys::System,n_m::Int64=4)
        @warn "Boundary conditions and load assumed to be 0"

        evals, evecs = (0,0)
        while any(evals .<= 0)
            evals, evecs = real.(Arpack.eigs(sys.Me,sys.Se,nev=n_m))
        end

        freqs = evals.^(-0.5)
        modes = u_to_Vh.(Ref(sys.problem.grid),eachcol(evecs))
        
        return evals, evecs, freqs, modes
    end

    function solve_dy_eigen(sys::System,n_m::Int64=4)
        evals, evecs, freqs, modes = get_vibrations(sys,n_m) 
        
        X(x::Float64) = [mode(x) for mode in modes]

        function get_T(IC::Matrix{Float64})
            @assert size(IC) == (sys.shape[2],2) "Wrong IC size for given system"

            as = evecs\IC[:,1]
            bs = (evecs\IC[:,2])./freqs

            T(t::Float64) = as.*cos.(freqs.*t).+bs.*sin.(freqs.*t)

            return T
        end

        return X, get_T
    end
end