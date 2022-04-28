using Flower1D, Plots, LaTeXStrings

function correctionFunctions()
    P = 5 # Order of the correction polynomials
    C = 2 # C scheme choice in vcjh

    x = collect(LinRange(-1.0, 1.0, 1000))
    h = vcjh(P, x, C = C)
    dh = vcjh(P, x; C = C, derivativeChoice = true)

    a = plot()
    plot(x, h[:, 1], label = L"h_L(x)", color = :black)
    plot!(x, h[:, 2], label = L"h_R(x)", color = :black, linestyle = :dash)

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
        xlims = (-1.0, 1)
    )
    savefig(string(@__DIR__) * "/h.pdf")

    a = plot()
    plot!(x, dh[:, 1], label = "\\partial h_L(x)" |> latexstring, color = :black)
    plot!(
        size = (600, 500),
        xlims = (-1.0, 1)
    )
    savefig(string(@__DIR__) * "/dhL.pdf")
end

correctionFunctions()