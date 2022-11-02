module BasisFunctions
    using Plots
    using LaTeXStrings

    ϕ₁(x) = 1 - 3x^2 + 2x^3
    ϕ₂(x) = x*(x - 1)^2
    ϕ₃(x) = 3x^2 - 2x^3
    ϕ₄(x) = x^2 *(x - 1)

    function form_functions()
        x = 0:0.01:1
        y1 = ϕ₁.(x)
        y2 = ϕ₂.(x)
        y3 = ϕ₃.(x)
        y4 = ϕ₄.(x)
        p = plot(title="Form functions", legendfontsize=15)
        plot!(p, x, y1, linewidth=3, color="red",label=L"\bar{\phi}_1")
        plot!(p, x, y2, linewidth=3, color="blue",label=L"\bar{\phi}_2")
        plot!(p, x, y3, linewidth=3, color="black",label=L"\bar{\phi}_3")
        plot!(p, x, y4, linewidth=3, color="orange",label=L"\bar{\phi}_4")
        savefig(p, "img/single/form_functions.svg")
    end

    function basis_functions()
        x = -0.3:0.001:0.3
        y1 = zeros(size(x))
        y2 = zeros(size(x))
        for (i,xi) in enumerate(x)
            if -0.1 < xi < 0
                y1[i] = ϕ₃((xi + 0.1)/0.1)
                y2[i] = 0.9*ϕ₄((xi + 0.1)/0.1)
            elseif 0 ≤ xi < 0.1
                y1[i] = ϕ₁(xi/0.1)
                y2[i] = 0.9*ϕ₂(xi/0.1)
            end
        end
        p = plot(title="Basis functions", legendfontsize=15)
        plot!(p, x, y1, linewidth=3, color="blue", label=L"\phi_{2i-1}")
        plot!(p, x, y2, linewidth=3, color="red", label=L"\phi_{2i}")
        savefig(p, "img/single/basis_functions.svg")
    end
end
