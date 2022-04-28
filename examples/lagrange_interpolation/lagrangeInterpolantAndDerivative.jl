using Revise, Flower1D
using Plots, LaTeXStrings

gr()

"""
Checks the correct implementation of the Lagrange interpolation as well as its derivative for
a given `f :: Function` and its derivative `df :: Function`.
"""
function lagrangeInterpolantAndDerivative(f :: Function, df :: Function)
    x = collect(LinRange(-1, 1, 100))
    plot(x, f.(x), label = "Target", color = :black)

    standardElement = StandardElement(1, 5, "legendre")
    fj = f.(standardElement.innerPoints)
    wj = lagrangeBarycentricWeights(standardElement.innerPoints)
    Pn = lagrangeInterpolationAt(x, standardElement.innerPoints, fj, wj)
    plot!(x, Pn, label = "Lagrange interpolant", linestyle = :dash)

    plot!(x, df.(x), label = "Target derivative")
    Dij = lagrangeDerivativeMatrix(standardElement.innerPoints)
    dPn = Dij * fj
    plot!(standardElement.innerPoints, dPn, label = "Lagrange interpolant derivative at nodes", seriestype = :scatter)

    # Format and save fig
    Plots.default(
        fontfamily = "Computer Modern",
        linewidth = 3,
        framestyle = :box,
        label = nothing,
        grid = false,
        margin=(Plots.Measures.Length(:mm, 2))
    )
    plot!(
        size = (600, 500),
        xlims = (-1, 1),
        minorticks=0.05,
        legend = :bottomright
    )
    savefig(string(@__DIR__) * "/lagrangeInterpolantAndDerivative.pdf")
end

f(x) = x ^ 2 + 3 * x + 1
df(x) = 2 * x + 3
f(x) = sin(x)
df(x) = cos(x)

lagrangeInterpolantAndDerivative(f, df)
