include("Beam1D.jl")
import Plots, LinearAlgebra
using LaTeXStrings

foeppl(xi, alpha, n)=  ifelse.(xi .> alpha, (xi .-alpha).^n, 0 )

function case_1_supported_beam_(l , q, EI)
	Xi(x) = x./l

	BCs  = Dict((0,'H')=>0,
                (0,'G')=>(q*l^2/24),
	            (1,'H')=>0,
                (1,'G')=>-(q*l^3/24))

	q_func(x) = q
    return x-> (q*Xi(x).^4 .- 2q*l*Xi(x).^3 .+ q*l^3*Xi(x)) / (24EI), BCs, q_func
end

function case_7_constant_load(l , q ,EI)
	q_func(x) = q
	Xi(x) = x./l 

	BCs  = Dict((0,'H')=>0,
	            (0,'G')=>0,
	            (1,'M')=>0,
	            (1,'Q')=>0)

	return  x-> (q*l^4/24 * (6*Xi(x).^2 - 4*Xi(x).^3 + Xi(x).^4)) / EI, BCs, q_func
end

function case_8_partly_constant_load(l, q, a,  EI)
	q_func(x) = ifelse.(x.>= a, q, 0)
	BCs  = Dict((0,'H')=>0,
	            (0,'G')=>0,
	            (1,'M')=>0,
	            (1,'Q')=>0)

	Xi(x) = x./l 
	alpha = a/l
	beta  = (l-a) /l
	return x-> q*l^4/24 * (foeppl(Xi(x), alpha, 4) - 4*beta*Xi(x).^3+ 6*beta*(2-beta)*Xi(x).^2) / EI, BCs, q_func
end

function case_9_decreasing_load(l , q ,EI)
    q_func(x) = q-q*x
	BCs  = Dict((0,'H')=>0,
	            (0,'G')=>0,
	            (1,'M')=>0,
	            (1,'Q')=>0)
	Xi(x) = x./l 
	return  x-> q*l^4/120 * (10*Xi(x).^2-10*Xi(x).^3 + 5*Xi(x).^4 -Xi(x).^5) / EI, BCs, q_func
 end

function case_10_free_momentum_at_a(l, M, a, EI)
	q_func(x) = 0
	BCs  = Dict((0,'H')=>0,
	            (0,'G')=>0,
	            (1,'M')=>M,
	            (1,'Q')=>0)
	Xi(x) = x./l 
	alpha = a/l
	return x-> -M*l^2/2 .* (Xi(x).^2 - foeppl(Xi(x),alpha, 2)) / EI, BCs, q_func
end

function constant_vs_decreasing_load()
    L = 1.0
    q = -3
    EI_const = 1
    grid = collect(LinRange(0,L,20))
    cases = Dict(
        "constant_load"  =>Dict("fnc"=>case_7_constant_load,  "input"=>(L,q,EI_const),"legend"=>:bottomleft,"color"=>"blue","dotcolor"=>"red"),
        "decreasing_load"=>Dict("fnc"=>case_9_decreasing_load,"input"=>(L,q,EI_const),"legend"=>:bottomleft,"color"=>"red","dotcolor"=>"blue"),
    )
    p = Plots.plot()
    for (key,value) in cases
        analytic_sol, BCs, q_func = value["fnc"](value["input"]...)
        # Parameters 
        pars = (mu=x->1 ,EI=x->EI_const, q=q_func)

        # Build and solve
        problem = Beam1D.Problem(pars, BCs, grid)
        sys = Beam1D.System(problem)


        u_numeric = sys.Se\sys.qe # 4 boundary conditions at end
        u_analytic = analytic_sol(grid)

        xs,ys = Beam1D.u_to_Vh(grid,u_numeric)

        # Plot
        fancy_key = uppercasefirst(replace(key,"_"=>" "))
        Plots.plot!(Beam1D.eval(xs,0.5),Beam1D.eval(ys,0.5),linewidth=3,label=fancy_key,color=value["color"])
        Plots.plot!(legend=value["legend"])
        Plots.ylabel!("w(x)")
        Plots.xlabel!("x")
    end
    Plots.savefig("img/single/constant_vs_decreasing_load.svg")
end

function test_all_cases()
    L = 1.0
    q = -3 # Negative y-direction
    M = -9
    a = 1
    EI_const = 1
    grid = collect(LinRange(0,L,20))

    cases = Dict(
        "supported_beam"        => Dict("fnc"=>case_1_supported_beam_,      "input"=>(L,q,  EI_const), "legend"=>:top),
        "constant_load"         => Dict("fnc"=>case_7_constant_load,        "input"=>(L,q,  EI_const), "legend"=>:bottomleft),
        "partly_constant_load"  => Dict("fnc"=>case_8_partly_constant_load, "input"=>(L,q,a,EI_const), "legend"=>:bottomleft),
        "decreasing_load"       => Dict("fnc"=>case_9_decreasing_load,      "input"=>(L,q,  EI_const), "legend"=>:bottomleft),
        "free_momentum_at_a"    => Dict("fnc"=>case_10_free_momentum_at_a,  "input"=>(L,M,a,EI_const), "legend"=>:topleft),
    )

    for (key,value) in cases
        analytic_sol, BCs, q_func = value["fnc"](value["input"]...)
        # Parameters 
        pars = (mu=x->1 ,EI=x->EI_const, q=q_func)

        # Build and solve
        problem = Beam1D.Problem(pars, BCs, grid)
        sys = Beam1D.System(problem)

        u_numeric = sys.Se\sys.qe # 4 boundary conditions at end
        u_analytic = analytic_sol(grid)

        word_print = uppercasefirst(replace(key,"_"=>" "))
        println(word_print * ", L2 norm: \n\t", LinearAlgebra.norm(u_numeric[1:2:end-4]- u_analytic))

        xs,ys = Beam1D.u_to_Vh(grid,u_numeric)

        # Plot
        points = [0,1]
        legend = hcat("Numerical solution",fill("",1,length(points)-1))
        p = Plots.plot()
        Plots.plot!(grid,u_analytic, label = "Analytical solution",color="blue")
        Plots.plot!(Beam1D.eval(xs,points),Beam1D.eval(ys,points),seriestype=:scatter,label=legend,
            color="red")
        Plots.plot!(legend=value["legend"])
        Plots.ylabel!("w(x)")
        Plots.xlabel!("x")
        Plots.savefig("img/single/"*key*".svg")
    end
end

function convergence_testing()
    L = 1.0
    q = -3
    EI_const = 1
    resolutions = Int.(floor.(10 .^(LinRange(log10(3),log10(1e3),50))))
    errors = zeros(length(resolutions))

    analytical_sol, BCs, q_func = case_1_supported_beam_(L,q, EI_const)
    pars = (mu=x->1 ,EI=x->EI_const, q=q_func)

    interpolation_points = 100
    hi = L / interpolation_points
    for (j,resolution) in enumerate(resolutions)
        h = L / resolution
        grid = collect(LinRange(0,L,resolution+1))

        problem = Beam1D.Problem(pars, BCs, grid)
        sys = Beam1D.System(problem)
        u_numeric = sys.Se\sys.qe
        u_analytic = analytical_sol(grid)

        xs,ys = Beam1D.u_to_Vh(grid, u_numeric)

        errors[j] = LinearAlgebra.norm(u_numeric[1:2:end-4]- u_analytic)
    end

    p = Plots.plot()
    Plots.plot!(resolutions, errors, xaxis=:log, yaxis=:log, seriestype=:scatter, 
        markersize=4, color="blue",label=false,
        title=L"Error wrt. $L_2$-norm with 100 interpolation points")
    Plots.ylabel!(L"\vert\vert{w_{\textrm{numerical}} - w_{\textrm{analytical}}}\vert\vert")
    Plots.xlabel!("# of elements")
    Plots.savefig("img/single/convergence.svg")
    Plots.savefig("report/overleaf/plots/convergence.svg")
    Plots.savefig("presentation/plots/convergence.svg")
end
