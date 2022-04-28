using Revise, Flower1D
using Plots, LaTeXStrings

gr()

function lagrangeMeshElements()
    f(x) = 1 / (1 + 16 * x^2)

    # Define constants
    T = Float64 # define precision
    N = 5 # number of elements per dimension
    E = 1 # element type, (1: line, 2: quad, 3: triangle)
    P = 2 # element order
    C = 1 # correction function parameter
    LMIN, LMAX = -1.0, 1.0 # domain boundary
    FTP = 100 # reference function test points
    ETP = 10 # element test points

    # Create mesh
    vertexPoints = collect(LinRange(LMIN, LMAX, N + 1))
    elementsVertices = Matrix{T}(undef, N, 2)
    for i ∈ 1:length(vertexPoints) - 1
        elementsVertices[i, :] = [vertexPoints[i], vertexPoints[i + 1]]
    end
    mesh = Mesh(E, P, "legendre", elementsVertices, C; T = T)

    # Plot target function
    x = collect(LinRange(LMIN, LMAX, FTP))
    plot(x, f.(x), label="f(x)", color="black")

    # Discrete element solution
    for (i, element) in enumerate(mesh.elements)
        xi = element.innerPoints # solution points
        yi = f.(xi) # solution values

        x = collect(LinRange(element.vertexPoints[1], element.vertexPoints[end], ETP))
        wi = lagrangeBarycentricWeights(xi)
        y = lagrangeInterpolationAt(x, xi, yi, wi)
        plot!(x, y, label="Element $i", markershape=:circle, linewidth = 2)
    end

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
        ylims = (0, 1.05),
        minorticks=0.05,
        legend = :topright
    )
    savefig(string(@__DIR__) * "/lagrangeMeshElements.pdf") # Saves the plot from p as a .pdf vector graphic
end

lagrangeMeshElements()