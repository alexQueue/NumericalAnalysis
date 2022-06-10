"""
Several beam bending examples for analytic solutions
"""

foeppl(xi, alpha, n)=  ifelse.(xi .> alpha, (xi .-alpha).^n, 0 )

function case_1_supported_beam_(l , a, F ,EI)
    alpha = a/l
    beta  = (l-a) /l
    Xi(x) = x./l 
    BCs  = Dict((0,'H')=>0,
                (0,'G')=>(((F*l^2)/6)*(beta-beta^3))/EI,
                (1,'H')=>0,
                (1,'G')=>-(((F*l^2)/6)*(alpha-alpha^3))/EI            
    )
    q_func(x) = q0#ifelse.(x.== a, q0, 0) #TODO we need F and not q
   return  x-> ((((F*l^3)/6).*(beta*Xi(x) .* (1-beta^2 .- Xi(x).^2) .+ foeppl(Xi(x), alpha, 3)))/EI), BCs, q_func 
end

function case_7_constant_load(l , q0 ,EI)
    q_func(x) = q0
    BCs  = Dict((0,'H')=>0,
    (0,'G')=>0,
    (1,'M')=>0,
    (1,'Q')=>0)

   Xi(x) = x./l 
   return  x-> (((q0*l^4)/24)* (6*Xi(x).^2 - 4*Xi(x).^3 + Xi(x).^4 ))/EI, BCs, q_func
end

function case_8_partly_constant_load(l, q0, a,  EI)
    q_func(x) = ifelse.(x.>= a, q0, 0)
    BCs  = Dict((0,'H')=>0,
    (0,'G')=>0,
    (1,'M')=>0,
    (1,'Q')=>0)

    Xi(x) = x./l 
    alpha = a/l
    beta  = (l-a) /l
    return x-> ((q0*l^4)/24)* (foeppl(Xi(x), alpha, 4) - 4*beta*Xi(x).^3+ 6*beta*(2-beta)*Xi(x).^2)/EI, BCs, q_func
end

function case_9_decreasing_load(l , q0 ,EI)
    q_func(x) = q0-q0*x
    BCs  = Dict((0,'H')=>0,
    (0,'G')=>0,
    (1,'M')=>0,
    (1,'Q')=>0)
    Xi(x) = x./l 
    return  x->  (q0*l^4)/120 * (10*Xi(x).^2-10*Xi(x).^3 + 5*Xi(x).^4 -Xi(x).^5)/EI, BCs, q_func
 end

 function case_10_free_momentum_at_a(l, a, M_0, EI)
    println("Free momentum set at point ", a)
    q_func(x) = 0
    BCs  = Dict((0,'H')=>0,
            (0,'G')=>0,
            (1,'M')=>M_0,
            (1,'Q')=>0)
    Xi(x) = x./l 
    alpha = a/l
    return x-> ((-M_0*l^2)/2 .* (Xi(x).^2- foeppl(Xi(x),alpha, 2)))/EI, BCs, q_func
 end
