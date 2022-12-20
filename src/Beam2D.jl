# Beam2D.jl
# Authors: Alex Quinlan, Pia Callmer, Nicola Sabbadini & Henry Jacobson
#
# TU Berlin
# Project Numerical Analysis SoSe 2022

module Beam2D
	using SparseArrays, LinearAlgebra, Plots #Stdlib imports
	import IterTools, Arpack, Printf, SciPy #External imports

    """
        struct Node

    Struct that holds information on which type of node, the coordinates as well (if relevant) the 
    force/moment/movable_direction

    The constructor doesnt populate the connecting_edges member since it has to be populated using
    edges construct, which in turn needs the specific node. (This isnt circular since an edge doesnt
    need to know which connecting edges a node has).
    """
    mutable struct Node
        type    ::String # FREE/FIXED/FORCE/MOVABLE
        coord   ::Vector{Float64} # 2 element vector of global cooridnates.
        number  ::Int64 # The iteration when this node object was created

        # The specific values for a node of type FORCE/MOVABLE
        force               ::Vector{Float64}
        moment              ::Float64
        movable_direction   ::Union{Vector{Float64},String}

        connecting_edges    ::Vector{Int64} # The edges that are connected to this node

        function Node(
                    type              ::Union{String,SubString{String}},
                    coord             ::Vector{Float64},
                    number            ::Int64;
                    force             ::Vector{Float64} = [0.0,0.0],
                    moment            ::Float64 = 0.0,
                    movable_direction ::Union{Vector{Float64},String} = "ALL"
                )
            if movable_direction != "ALL"
                movable_direction = movable_direction / norm(movable_direction) # Normalize the movable direction
            end
            new(type,coord,number,force,moment,movable_direction,[])
        end
    end

    """
        struct Edge

    Struct that holds information on which nodes an edge is connected to and the
    relevant finite element method information.
    """
    mutable struct Edge
        nodes       ::Vector{Node} # 2 element vector which nodes the edge is construed from
        grid        ::Vector{Float64} # The finite element gridpoints used in local coordinates
        gridlen     ::Int64 # How many elements are in the grid
        len         ::Float64 # The length of the edge
        index_start ::Int64 # Which index in the whole framework the edge's start is

        # Functions of x and t for the parameters
        E   ::Function
        I   ::Function
        A   ::Function
        mu  ::Function
    end

    """
        struct Problem

    Struct that holds information of the framework of nodes and edges 
    """
    struct Problem
        nodes ::Vector{Node} # Collection of the nodes in the framework
		edges ::Vector{Edge} # Collection of the edges in the framework
        size  ::Int64 # Amount of elements in the whole framework

        function Problem(file::String)
            nodes, edges = problem_constructor(file)

            last_edge = edges[end]
            size = last_edge.index_start + length(last_edge.grid)*3 - 1

            # Add the corresponding edges for each node
            for (i, edge) in enumerate(edges)
                for node in edge.nodes
                    push!(node.connecting_edges, i)
                end
            end

            new(nodes,edges,size)
        end
	end

    """
        problem_constructor(file::String)

    Constructs an instance of the Problem struct from `file`.
    """
    function problem_constructor(file::String)
        # Read file and save initial lists
        io = open(file, "r")
        txt = read(io, String)
        fl_as_list = split(txt,"\n")
        fl_as_list = strip.(fl_as_list) # Remove trailing whitespace
        filter!(x->x!= "", fl_as_list) # Remove empty elements
        filter!(x->x[1] != '#', fl_as_list) # Remove line starting with a comment
        fl_as_list = String.(fl_as_list) # Cast from SubString to String
        fl_as_list = striplinecomment.(fl_as_list) # Remove comments
        indices = findall(x->x in ["NODES","EDGES","TYPE","PARAMETERS"], fl_as_list) 
        splitted = getindex.(
             Ref(fl_as_list), 
             UnitRange.(
                [1; indices .+ 1],
                [indices .- 1; length(fl_as_list)]
            )
        ) # Split into 4 lists
        nodes,edges,types,params = splitted[2:end] # 2:end because 1st element is empty as we split on "NODES"

        node_data = split.(nodes, " ")
        nodes = [parse.(Float64, x) for x in node_data] # Get the coordinates of the nodes to the nodes list

        edges = edges_setup(edges)
        
        types = get_node_type_list(types)

        # Create and populate a list of Node objects
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

        parameters = Dict()
        for param in params
            str,val = split(param, " ", limit=2)
            f = eval(Meta.parse("x -> " * val)) # Evaluates the string to julia code
            parameters[str] = x -> Base.invokelatest(f,x)
        end
        E = parameters["E"]; I = parameters["I"]; A = parameters["A"]; mu = parameters["mu"]

        # Create and populate a list of Edge objects
        index_cnt = 1
        Edges = Vector{Edge}(undef, length(edges))
        for (i,edge) in enumerate(edges)
            L = norm(Nodes[edge[1]].coord - Nodes[edge[2]].coord)
            gridpoints = edge[3]
            grid = collect(LinRange(0.0,L,gridpoints)) # 0 to length of beam

            if edge[4] == []
                params = [E,I,A,mu]
            else
                params = edge[4]
            end

            Edges[i] = Edge([ Nodes[edge[1]], Nodes[edge[2]] ], grid, length(grid), L, index_cnt, params...)
            index_cnt += length(grid)*3 # v_1 -> v_n, and w_1 -> w_2n makes 3n
        end

        return Nodes,Edges
    end
    
    """
        get_node_type_list(types::Vector)

    Returns the force+moment or movable direction, or just nothing if free/force from the
    `types` vector.
    """
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

    """
        edges_setup(edges_data::Vector)

    Parses each element in `edges_data` to return a vector of information
    on nodes, # of gridpointts and parameter expressions.
    """
    function edges_setup(edges_data::Vector)
        data = split.(edges_data, " ", limit=3) # Create a list of each edge line, group gp and params together
        edges = []

        # Regex for finding specific parameters / gridpoints
        re_params = r"params=\[(.*),(.*),(.*),(.*)\]"
        re_gp = r"gp=(\d*)"

        for edge_data in data
            nodes = [parse(Int, x) for x in edge_data[1:2]]
            if length(edge_data) == 2
                push!(edge_data,"")
            end
            # Find the amount of gridpoints / parameter expressions specified (if at all)
            gp_match = match(re_gp, edge_data[3])
            if gp_match !== nothing
                gridpoints = eval(Meta.parse(gp_match.captures[1]))
            else
                gridpoints = 2
            end
            params_match = match(re_params, edge_data[3])
            if params_match !== nothing
                params = []
                # params is a list of each specific parameter function
                for capture in params_match.captures
                    f = eval(Meta.parse("x -> " * capture))
                    param = x -> Base.invokelatest(f,x)
                    push!(params, param)
                end
            else
                params = []
            end
            push!(edges, [nodes...,gridpoints,params])
        end
        return edges
    end

    # Copied from internet obviously
    """
        striplinecomment(a::String, cchars::String="#")

    Strips comment `cchars` from string `a`.
    Copied from the internet of course.
    """
    function striplinecomment(a::String, cchars::String="#")
        b = strip(a)
        0 < length(cchars) || return b
        for c in cchars
            r = Regex(Printf.@sprintf "\\%c.*" c)
            b = replace(b, r => "")
        end
        strip(b)
    end

    """
        constraints_construction(problem::Problem)

    Constructs the constraint matrix C and constraint forces f from `problem`
    object of class Problem.
    """
    function constraints_construction(problem::Problem)
        # Count size of cᵀ, denoted r
        r = 0
        for node in problem.nodes
            n_cnct_edges = length(node.connecting_edges)
            if node.type == "FIXED"
                r += n_cnct_edges*3 # 3 conditions per connecting edge
            elseif node.type == "MOVABLE"
                r += 1 # one bearing condition for edge 1
                if n_cnct_edges >= 2
                    r += (n_cnct_edges - 1)*3
                end
            else
                r += (n_cnct_edges - 1)*3 # 3 conditions per pair of edges
            end
        end

        # Initialize sparse empty matrices
        C = spzeros(problem.size,r)
        f = spzeros(problem.size)

        # Loop through all nodes and add constraints from each
        i = 1
        for node in problem.nodes
            if node.type == "FIXED"
                # Loop through the edges of the node and make each edge fixed
                # in the side of the node.
                for edge in problem.edges[node.connecting_edges]
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
                end
            elseif node.type == "MOVABLE"
                edge = problem.edges[node.connecting_edges[1]]
                phi = edge_angle(edge)

                # No movement in movable direction means the dot product of movement
                # vector with rotated displacement vector 90^o is 0
                j = linking_index(edge, node)
                dx = node.movable_direction[1]; dy = node.movable_direction[2]
                
                rot = [cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)]
                m_orth = (rot*[dx;dy])'
                v_mov = m_orth*[cos(phi);sin(phi)]
                w_mov = m_orth*[-sin(phi);cos(phi)]

                C[j,i] = [v_mov, w_mov]
                i += 1

                # Even movable nodes can have more than one connecting edge 
                # and then we need stiffness and linking again
                i = connecting_edges_conditions!(problem, node, C, i)
            else # "FORCE/FREE"
                if node.type == "FORCE"
                    edge = problem.edges[node.connecting_edges[1]]
                    angle = -edge_angle(edge) # Rotation is done counter-clockwise to the angle
                    j1,j2 = linking_index(edge, node)
                    fx = node.force[1]; fy = node.force[2]

                    # Add force in correct direction on the local coordinate system of node
                    f[[j1,j2]] = [cos(angle) -sin(angle); sin(angle) cos(angle)]*[fx; fy] 
                end
                i = connecting_edges_conditions!(problem, node, C, i)
            end
        end
        C,f
    end

    """
        connecting_edges_conditions!(problem::Problem, node::Node, C, i::Int64)

    Inserts the linking and stiffness condition for the node in the constraint matrix,
    returns the iteration i after its done.
    """
    function connecting_edges_conditions!(problem::Problem, node::Node, C, i::Int64)
        first_edge = problem.edges[node.connecting_edges[1]] # Compare everything to first edge (arbitrary which)
        phi_1 = edge_angle(first_edge)
        
        for rem_edge in problem.edges[node.connecting_edges[2:end]]
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

    """
        edge_node_order(edge, node)

    Return "first" if `edge` starts at `node`.
    """
    function edge_node_order(edge, node)
        edge.nodes[1].number == node.number ? "first" : "last"
    end

    """
        fixd_index(edge, node)

    Returns thes indices of the constraints of a FIXED node.
    """
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
                edge.index_start + edge.gridlen*3-2,
                edge.index_start + edge.gridlen*3-1
            ]
        return inds
    end

    """
        linking_index(edge, node)

    Returns the indices for the linking condition for the `edge` with order "last" or "first".
    It is the same indices for the force.
    """
    function linking_index(edge, node)
        order = edge_node_order(edge, node)
        inds = order == "first" ? 
            [edge.index_start, 
                edge.index_start + edge.gridlen] :
            [edge.index_start + edge.gridlen - 1,
                edge.index_start + edge.gridlen*3 - 2]
        return inds
    end

    """
        stiffness_index(edge, node)

    Returns the index corresponding to a stiff node.
    """
    function stiffness_index(edge, node)
        order = edge_node_order(edge, node)
        index = order == "last" ? 
            edge.index_start + edge.gridlen*3 - 1 :
            edge.index_start + edge.gridlen + 1 
        return index
    end

    """
        edge_angle(edge::Edge)

    Returns the angle of an edge to the x-axis.
    """
    function edge_angle(edge::Edge)
        dy = edge.nodes[2].coord[2] - edge.nodes[1].coord[2]
        dx = edge.nodes[2].coord[1] - edge.nodes[1].coord[1]
        atan(dy,dx)
    end

    """
        struct System

    Struct that contains the matrices that define the problem, as well
    as the shape which is the size of the S/M matrix and the length of the force vector 
    (the same as the # of columns of the contraint matrix).
    """
	struct System
		problem ::Problem
        shape   ::Tuple{Int,Int}
		Me      ::SparseMatrixCSC{Float64,Int64}
		Se      ::SparseMatrixCSC{Float64,Int64}
		qe      ::Vector{Float64}

        """
            System(problem::Problem)

        Constructor of a `System` object from a `problem`.
        """
		function System(problem::Problem)
            n = sum([edge.gridlen*3 for edge in problem.edges])

            C,f = constraints_construction(problem)
            r = size(C)[2]

            Me = spzeros(n+r,n+r)
            Se = spzeros(n+r,n+r)
            qe = zeros(Float64,n+r)

            Se[1:n,n+1:n+r] = C
            Se[n+1:n+r,1:n] = Transpose(C)
            qe[1:n] = f

            # Physics
            # Shape functions and derivatives on element [0,h], with t in [0,1]
            phi_T_0(h,t) = [ t^2*(2*t-3)+1, t*(t-1)^2*h  ,-t^2*(2*t-3)  , t^2*(t-1)*h   ]
            phi_T_1(h,t) = [ 6*t*(t-1)/h  , (t-1)*(3*t-1),-6*t*(t-1)/h  , t*(3*t-2)     ] # Not used
            phi_T_2(h,t) = [ (12*t-6)/h^2 , (6*t-4)/h    ,-(12*t-6)/h^2 , (6*t-2)/h     ]
            phi_T_3(h,t) = [ 12/h^3       , 6/h^2        ,-12/h^3       , 6/h^2         ] # Not used

            phi_L_0(h,t) = [ t            , 1-t]
            phi_L_1(h,t) = [ 1/h          ,-1/h]

            # 1D, 3 point Gaussian quadrature on element [0,h], with t in [0,1]
			GQ3(h,f) = 5*h/18*f(1/2-sqrt(3/20)) +
                       4*h/9 *f(1/2           ) +
                       5*h/18*f(1/2+sqrt(3/20))

            # Longitudinal - Transversal ordering
            M_base_fncs = [phi_L_0, phi_T_0]
            S_base_fncs = [phi_L_1, phi_T_2]

            # Local System for S for element [o,o+h], with t in [0,1]
            S_loc(h,o,base_fnc,E_) = GQ3(h,t -> base_fnc(h,t)*base_fnc(h,t)'*E_(o+h*t))
            
            for edge in problem.edges
                EI(x) = edge.E(x) * edge.I(x)
                EA(x) = edge.E(x) * edge.A(x)
                E_s = [EA,EI]

                # Local System for M for element [o,o+h], with t in [0,1]
                M_loc(h,o,base_fnc) = GQ3(h,t -> base_fnc(h,t)*base_fnc(h,t)'*edge.mu(o+h*t))

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
                    
                    # Loop through the local indices
                    for (i,h,o) in zip(collect.(IterTools.partition(1:n,2j,j)), diff(edge.grid), edge.grid)
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
        plot_static(sys::System)

    Solves the problem defined in `sys` and returns a plot object of the plot.
    """
    function plot_static(sys::System)
        u = sys.Se\sys.qe
        xs,ys = u_to_Vh(sys.problem, u)
        xs_undeformed,ys_undeformed = u_to_Vh(sys.problem, zeros(size(u)...))
        p = plot()
        plot!(xs_undeformed,ys_undeformed,0,1,color="black",label=false,linewidth=2,linestyle=:dot)
        plot!(xs,ys,0,1,color="black",label=false,linewidth=2)
        p
    end
    
    function animate_framework(sys::System, savelocation, tf=100, tsteps=100)
        get_sol = Beam2D.solve_dy_eigen(sys)
        u = sys.Se\sys.qe

        IC = [u zeros(size(u)...,2)]
        sol = get_sol(IC)
        ts = collect(LinRange(0,tf,tsteps,))

        undeformed = Beam2D.u_to_Vh(sys.problem,zeros(size(u)...))

        anim = Plots.@animate for t in ts
            Plots.plot(
                undeformed[1],undeformed[2],0,1,
                color="black",label=false,linewidth=2,linestyle=:dot
            )
            Plots.plot!(sol(t)...,0,1,color="black",label=false,linewidth=2)
        end
        Plots.gif(anim, savelocation, fps=15)
    end

    """
        u_to_Vh(problem::Problem,u::AbstractVector{Float64})

    Converts a solution `u` to parametric curves using `problem` object.
    Returns two vectors of the x-functions and the y-functions, parametrized
    from 0 to 1.
    """
    function u_to_Vh(problem::Problem, u::AbstractVector{Float64}) 
        phi(t) = [
            t^2*(2*t-3)+1,
            t*(t-1)^2,
            -t^2*(2*t-3),
            t^2*(t-1)
         ] # The basis functions

        # Number of parametric curves, one for each element on each edge
        n_elements = sum([edge.gridlen-1 for edge in problem.edges]) 

        xs = Vector{Function}(undef,n_elements)
        ys = Vector{Function}(undef,n_elements)

        i = 1
        for edge in problem.edges
            edge_dir = edge.nodes[2].coord - edge.nodes[1].coord # Direction of the edge in global coordinates
            edge_start = edge.nodes[1].coord
            pos(t) = edge_start + t*edge_dir # Function to get the position of the start finite element on the edge

            rot_matrix = [edge_dir[1] -edge_dir[2]; edge_dir[2] edge_dir[1]]./norm(edge_dir)
            h = edge.len / edge.gridlen # Constant length elements for each edge

            for element_nr in 1:edge.gridlen-1
                # Get the positions of the start and end of the element in the global coordinate system for adding to the deformed values
                pos1 = pos((element_nr-1)/(edge.gridlen-1))
                pos2 = pos(element_nr/(edge.gridlen-1))

                # Get the indices in the global matrix for the specific elements
                v_inds = edge.index_start + element_nr - 1 .+ (0:1)
                w_inds = (edge.index_start + edge.gridlen) .+ (2*(element_nr-1):2*(element_nr-1)+3)
                (x0,x1,y0,g0,y1,g1) = u[[v_inds...,w_inds...]]

                q0_x,q0_y = pos1 .+ rot_matrix * [x0,y0]
                q1_x,q1_y = pos2 .+ rot_matrix * [x1,y1]
                u0_x,u0_y = rot_matrix * [1,g0] .* h
                u1_x,u1_y = rot_matrix * [1,g1] .* h

                # Determine the x and y values using Hermite polynomials
                xs[i] = t -> dot([q0_x,u0_x,q1_x,u1_x],phi(t))
                ys[i] = t -> dot([q0_y,u0_y,q1_y,u1_y],phi(t))
                i += 1
            end
        end

        return xs, ys
    end

    """
        solve_dy_Newmark(sys::System,IC::Matrix{Float64},times::Vector{Float64})

    Solves the problem using Newmarks method using the initial conditions `IC` for the `times` specified
    and returns a vector of the parametric curves for the whole framework for each timepoint.
    """
    function solve_dy_Newmark(sys::System, IC::Matrix{Float64}, times::Vector{Float64})
        N = sum(sys.shape)
        @assert size(IC) == (N,3) "Wrong IC size for given system"
        @assert length(times) >= 2 "Must have an initial and final time"
        @assert all(diff(times) .> 0) "Times must be ascending"
        
        n_t = length(times)
        
        u_0 = Array{Float64,2}(undef,N,n_t)
        u_1 = Array{Float64,2}(undef,N,n_t)
        u_2 = Array{Float64,2}(undef,N,n_t)
        
        # Put the initial condition in the first column of the solution
        u_0[:,1], u_1[:,1], u_2[:,1] = eachcol(IC)

        beta, gamma = (1/4,1/2) # Optimal parameters for the problem.

        for (t,h) in enumerate(diff(times))
            u_0s = u_0[:,t] + h*u_1[:,t] + (0.5-beta)*h^2*u_2[:,t]
            u_1s = u_1[:,t] + (1-gamma)*h*u_2[:,t]

            u_2[:,t+1] = (sys.Me+beta*h^2*sys.Se)\(sys.qe-sys.Se*u_0s)
            u_1[:,t+1] = u_1s + gamma*h*u_2[:,t+1]
            u_0[:,t+1] = u_0s + beta*h^2*u_2[:,t+1]
        end
        
        return [u_to_Vh(sys.problem,u) for u in eachcol(u_0)]
    end

    """
        get_vibrations(sys::System,n_m::Int64=4)

    Returns the first `n_m` generalized eigenvalues, eigenvectors, frequencies 
    and modes of the system matrices.
    """
    function get_vibrations(sys::System, n_m::Int64=4)
        @warn "Boundary conditions and load assumed to be 0"

        # Arpack generalized eigenvalues in julia didnt work for the 2D case for some reason
        # So we had to fall back to SciPy's implementation.
        Me = SciPy.sparse.csc_matrix(sys.Me)
        Se = SciPy.sparse.csc_matrix(sys.Se) 

        evals, evecs = real.(SciPy.sparse.linalg.eigs(Me,k=n_m,M=Se))

        freqs = evals.^(-0.5)
        modes = u_to_Vh.(Ref(sys.problem),eachcol(evecs))
        
        return evals, evecs, freqs, modes
    end

    """
        solve_dy_eigen(sys::System,n_m::Int64=4)

    Returns a function that returns a function for the eigenvalue method solution
    that takes in a given initial condition. Takes in a `sys` object of struct `System`
    and the number of eigenvaluees `n_m` to use.
    """
    function solve_dy_eigen(sys::System, n_m::Int64=4)
        evals, evecs, freqs, modes = get_vibrations(sys,n_m) 
        
        function get_sol(IC::Matrix{Float64})
            @assert size(IC) == (sum(sys.shape),3) "Wrong IC size for given system"

            as = evecs\IC[:,1]
            bs = (evecs\IC[:,2])./freqs

            sol(t) = u_to_Vh(sys.problem,evecs*(as.*cos.(freqs.*t).+bs.*sin.(freqs.*t)))

            return sol
        end

        return get_sol
    end
end
