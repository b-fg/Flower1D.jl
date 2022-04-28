using Revise, Flower1D
using Plots, LaTeXStrings

gr()

function checkLagrangeBasis()
    P = 3 # order
    NTP = 100 # number of test points
    standardElement = StandardElement(1, P, "legendre")
    out("standardElement", standardElement)

    xi = collect(LinRange(-1.0, 1.0, NTP))
    lj = zeros(eltype(xi), (P + 1, length(xi)))
    dlj = zeros(eltype(xi), (P + 1, length(xi)))
    for (i, x) ∈ enumerate(xi)
        lj[:, i] = lagrangeBasisAt(x, standardElement.innerPoints)
        dlj[:, i] = dLagrangeBasisAt(x, standardElement.innerPoints)
    end
    a = plot()
    for p ∈ 1:P + 1
        plot!(xi, lj[p, :], label = ("l_" * string(p) * "(x)") |> latexstring, color = palette(:tab10)[p])
        plot!(xi, dlj[p, :], label = ("\\partial l_" * string(p) * "(x)") |> latexstring, linestyle = :dash, color = palette(:tab10)[p])
    end

    @assert sum(lj) / NTP ≈ 1.0 "Error: the summation of the lagrange basis polynomials is not 1.0 at each location."

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
        legend = :topright
    )
    savefig(string(@__DIR__) * "/checkLagrangeBasis.pdf")
end

checkLagrangeBasis()