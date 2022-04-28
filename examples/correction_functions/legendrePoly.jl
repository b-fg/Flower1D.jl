using Flower1D, Plots, LaTeXStrings

function legendrePoly()
    P = 4 # Max order of the polynomials
    x = collect(LinRange(-1.0, 1.0, 1000))
    p = collect(1:P)
    Lp = pLegendreAt.(p, x')
    a = plot()
    for p ∈ 1:P
        plot!(x, Lp[p, :], label = ("L" * string(p)) |> latexstring, color = palette(:tab10)[p])
    end
    plot!(xlims = (-1, 1), ylims = (-6, 6))

    # Plot derivative
    dLp = dLegendreAt.(p, x')
    for p ∈ 1:P
        plot!(x, dLp[p, :], label = ("\\partial L" * string(p)) |> latexstring, color = palette(:tab10)[p], linestyle = :dash)
    end

    # Format and save fig
    Plots.default(
        fontfamily = "Computer Modern",
        linewidth = 2,
        framestyle = :box,
        label = nothing,
        grid = false,
        margin=(Plots.Measures.Length(:mm, 2))
    )
    plot!(
        size = (600, 500),
        xlims = (-1, 1),
        legend = :bottomright
    )
    savefig(string(@__DIR__) * "/legendrePoly.pdf")

end

legendrePoly()