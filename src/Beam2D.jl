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
            Edges[i] = Edge([ Nodes[edge[1]], Nodes[edge[2]] ], grid, index_cnt)
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

    # TODO - Add values to vector of values of constraints herein as well
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
                # v
                j = Problem.edges[node.connecting_edges[1]].index_start
                C[j,i] = 1
                i += 1
                
                # w
                j += length(Problem.edges[node.connecting_edges[1]].grid)
                C[j,i] = 1
                i += 1
                
                # w'
                j += 1
                C[j,i] = 1
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
                    order = edge_node_order(edge, node)
                    j1,j2 = linking_index(edge, order)
                    f[[j1,j2]] = [node.force[1]*cos(angle) ,node.force[2]*sin(angle)]
                end
                i = connecting_edges_conditions!(Problem, node, C, i)
            end
        end
        C,f
    end

    function edge_node_order(edge, node)
        edge.nodes[1].number == node.number ? "last" : "first"
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
            C[[j1,j2],i] = [1,-1]
            i += 1
            
            # Linking condition
            phi_2 = edge_angle(rem_edge)
            j1 = linking_index(first_edge, first_order)
            j2 = linking_index(rem_edge, rem_order)
            
            C[[j1...,j2...],i] = [cos(phi_1), -sin(phi_1), -cos(phi_2), sin(phi_2)]
            i += 1
            
            C[[j1...,j2...],i] = [sin(phi_1), cos(phi_1), -sin(phi_2), -cos(phi_2)]
            i += 1
        end
        return i
    end

    # Returns the indices for the linking condition for the edge with order "last" or "first"
    # Is the same indices for the force (i think)
    function linking_index(edge, order)
        idx = order == "first" ? 
            [edge.index_start, 
                edge.index_start + length(edge.grid)] :
            [edge.index_start + length(edge.grid) - 1,
                edge.index_start + length(edge.grid)*3 - 2]
        return idx
    end

    function stiffness_index(edge, order)
        idx = order == "first" ? 
            edge.index_start + length(edge.grid)*3 - 1 :
            edge.index_start + length(edge.grid) + 1 
        return idx
    end

    function edge_angle(edge::Edge)
        dy = edge.nodes[2].coord[2] - edge.nodes[1].coord[2]
        dx = edge.nodes[2].coord[1] - edge.nodes[1].coord[1]
        atan(dy,dx)
    end
    
    function edge_length(edge::Edge)
        dy = edge.nodes[2].coord[2] - edge.nodes[1].coord[2]
        dx = edge.nodes[2].coord[1] - edge.nodes[1].coord[1]
        return 
        
    end

	struct System
		problem ::Problem
        shape   ::Tuple{Int,Int}
		Me      ::SparseMatrixCSC{Float64,Int64}
		Se      ::SparseMatrixCSC{Float64,Int64}
		qe      ::Vector{Float64}

		function System(problem::Problem)
            n = sum([length(edge.grid)*3 for edge in problem.edges])

            C,f = C_matrix_construction(problem)
            r = size(C)[2]

            Me = spzeros(n+r,n+r)
            Se = spzeros(n+r,n+r)
            qe = zeros(Float64,n+r)

            Se[1:n,n+1:n+r] = C
            Se[n+1:n+r,1:n] = Transpose(C)
            qe[1:n] = f

            # Physics
            # Shape functions and derivatives on element [0,h], with p in [0,1]

            # Transversal eq
            phi_0_T(h,p) = [ p^2*(2*p-3)+1, p*(p-1)^2*h  ,-p^2*(2*p-3)  , p^2*(p-1)*h   ]
            phi_1_T(h,p) = [ 6*p*(p-1)/h  , (p-1)*(3*p-1),-6*p*(p-1)/h  , p*(3*p-2)     ]
            phi_2_T(h,p) = [ (12*p-6)/h^2 , (6*p-4)/h    ,-(12*p-6)/h^2 , (6*p-2)/h     ]
            phi_3_T(h,p) = [ 12/h^3       , 6/h^2        ,-12/h^3       , 6/h^2         ]

            # Longitudinal eq
            phi_0_L(h,p) = [p, 1 - p]
            phi_1_L(h,p) = [1/h, -1/h]

            # 1D, 3 point Gaussian quadrature on element [0,h], with p in [0,1]
			GQ3(h,f) = 5*h/18*f(1/2-sqrt(3/20)) +
                       4*h/9 *f(1/2           ) +
                       5*h/18*f(1/2+sqrt(3/20))

            EI(x) = problem.E(x) * problem.I(x)
            EA(x) = problem.E(x) * problem.A(x)

            # Longitudinal - Transversal ordering
            E_s = [EA,EI]
            M_base_fncs = [phi_0_L, phi_0_T]
            S_base_fncs = [phi_1_L, phi_2_T]

			# Local System for element [o,o+h], with p in [0,1]
			M_loc(h,o,base_fnc) = GQ3(h,p -> base_fnc(h,p)*base_fnc(h,p)'*problem.mu(o+h*p))
            S_loc(h,o,base_fnc,E_) = GQ3(h,p -> base_fnc(h,p)*base_fnc(h,p)'*E_(o+h*p))

            for edge in problem.edges
                partitions = [
                    # Longitudinal equation
                    edge.index_start:edge.index_start+length(edge.grid)-1,
                    # Transversal equation
                    edge.index_start+length(edge.grid):edge.index_start+length(edge.grid)*3-1,
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

    """
    Returns a vector with 4 points, 1 and 3 at the start and end of the two beams
    and 2 and 4 makes up the gradient
    """
    function to_global(problem::Problem, u::Vector{Float64})
        n = Int(sum([length(edge.grid)/2 for edge in problem.edges]))
        pos_vector = [[] for _=1:n]
        i = 1
        for edge in problem.edges
            OG_pos = edge.nodes[1].coord
            phi = edge_angle(edge)
            Dphi = [cos(phi) -sin(phi); sin(phi) cos(phi)]
            L = norm(edge.nodes[1].coord - edge.nodes[2].coord)
            

            Hermite2Bezier = [1 0 0 0; 0 0 0 1; -3 3 0 0; 0 0 -3 3]^-1

            m = Int(length(edge.grid)/2)
            for j in 1:m
                v_start = edge.index_start
                w_start = edge.index_start + length(edge.grid)
                idx_jump = (j-1)*6
                p1 = Dphi * [u[v_start + idx_jump]; u[w_start + idx_jump]]
                r1 = Dphi * [u[w_start + idx_jump + 1]; 0]

                p2 = Dphi * [u[v_start + idx_jump + 1]; u[w_start + idx_jump + 2]]
                r2 = Dphi * [u[w_start + idx_jump + 3]; 0]

                hermite_points = [p1 p2 r1 r2]'
                bezier_points = Hermite2Bezier * hermite_points

                pos_global = [x + OG_pos for x in eachrow(bezier_points)]
                pos_vector[i] = pos_global
                i += 1
            end
        end
        pos_vector
    end
end