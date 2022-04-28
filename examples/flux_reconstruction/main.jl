using Revise, Flower1D, Plots, LaTeXStrings

# Initial conditions
f1(x) = 1 / (1 + 16 * x ^ 2)
f2(x) = sin(2*pi * x)
f3(x) = exp(2.0 * x ^ 2)
f4(x) = sin(x)

# Constants
T = Float64 # define precision

PDE = 2 # Burgers flux: f = 0.5 * u * u
BC = 0 # periodic
LMIN, LMAX = -1, 1 # domain boundaries
N = 5 # number of elements
E = 1 # elements type
P = 2 # element order
C = 0 # correction function parameter
IC = Function[f1]
NP = 100 # number of points to plot the interpolant at each element

function main()
    # Create mesh
    vertexPoints = collect(LinRange(LMIN, LMAX, N + 1))
    elementsVertices = Matrix{T}(undef, N, 2)
    for i ∈ 1:length(vertexPoints) - 1
        elementsVertices[i, :] = [vertexPoints[i], vertexPoints[i + 1]]
    end
    mesh = Mesh(E, P, "legendre", elementsVertices, C; T = T)
    innerPoints = vcat([collect(transpose(mesh.elements[n].innerPoints)) for n ∈ 1:meshSize(mesh)]...)

    # Initialize flow
    flow = Flower1D.Flow(N, P, PDE, BC)
    init!(flow, IC, innerPoints)
    solutionBoundaryConditions!(flow)

    # Compute advection flux at solution points: flow.fδD
    @inbounds for m ∈ 1:size(flow.uδD, 2), n ∈ 1:size(flow.uδD, 1)
        flow.fδD[n, m, :] .= flux(flow.wδD[n, m, :], flow.PDE)
    end

    # Compute numerical fluxes at element interface: flow.fδI
    basisFace = mesh.standardElement.basisFace
    basisInner = mesh.standardElement.basisInner
    innerPointsξ = mesh.standardElement.innerPoints

    reconstructSolutionAtFace!(flow, basisFace)
    reconstructFluxesAtFace!(flow, basisFace)
    riemannFluxes!(flow)

    # Plot data
    p1 = plot()
    ξ = collect(LinRange(-1, 1, NP))
    for n ∈ 2:size(flow.uδD, 1) - 1
        x = xFromξ.(ξ, mesh.elements[n].facePoints ...) # element inner points in physical coords
        faceIndexL, faceIndexR = element2faces(n)

        # Plot reference function
        plot!(x, IC[1].(x), label = L"u(x)", color = :black)

        # Plot flow.uδD nodal values
        plot!(mesh.elements[n].innerPoints, flow.uδD[n, :, 1], color = :orange, seriestype = :scatter, label = "")

        # Plot flow.uδD(x) interpolant
        uδD = zeros(NP)
        for i ∈ 1:NP
            basisInner = lagrangeBasisAt(ξ[i], innerPointsξ)
            uδD[i] = sum(flow.uδD[n, :] .* basisInner)
        end
        plot!(x, uδD, label = "u^{\\delta D}(x)" |> latexstring, color = :orange)

        # Plot flow.fδD nodal values
        plot!(mesh.elements[n].innerPoints, flow.fδD[n, :, 1], color = :blue, seriestype = :scatter, label = "")

        # Plot flow.fδD(x) interpolant
        fδD = zeros(NP)
        for i ∈ 1:NP
            basisInner = lagrangeBasisAt(ξ[i], innerPointsξ)
            fδD[i] = sum(flow.fδD[n, :] .* basisInner)
        end
        plot!(x, fδD, label = "f^{\\delta D}(x)" |> latexstring, color = :blue)

        # Plot flow.fK (extrapolated flux values at faces) at element faces
        plot!([mesh.elements[n].facePoints[1]], [flow.fK[faceIndexL, 2]], label = "f^{K}" |> latexstring,
            seriestype = :scatter, color = :cyan)
        plot!([mesh.elements[n].facePoints[2]], [flow.fK[faceIndexR, 1]], label = "f^{K}" |> latexstring,
            seriestype = :scatter, color = :cyan)

        # Plot flow.fδI (common face flux obtained with a Riemann solver) at element faces
        plot!([mesh.elements[n].facePoints[1]], [flow.fδI[faceIndexL, 1]], label = "f^{\\delta I}" |> latexstring,
            seriestype = :scatter, color = :yellow)
        if n == N + 1
            plot!([mesh.elements[n].facePoints[2]], [flow.fδI[faceIndexR, 1]], label = "f^{\\delta I}" |> latexstring,
                seriestype = :scatter, color = :yellow)
        end

        # Plot correction flux fδC and continuos flux fδ
        fδC = zeros(NP)
        fδ = zeros(NP)
        for i ∈ 1:NP
            hL, hR = vcjh(P, ξ[i]; C = C)
            jumpL = flow.fδI[faceIndexL, 1] - flow.fK[faceIndexL, 2, 1]
            jumpR = flow.fδI[faceIndexR, 1] - flow.fK[faceIndexR, 1, 1]
            fδC[i] = jumpL * hL + jumpR * hR
            fδ[i] = fδC[i] + fδD[i]
        end
        plot!(x, fδC, label = "f^{\\delta C}(x)" |> latexstring, color = :green)
        plot!(x, fδ, label = "f^{\\delta}(x)=f^{\\delta D}+f^{\\delta C}" |> latexstring, color = :red, linestyle = :dash)
    end

    # Format and save fig
    # remove repeated labels in legend
    known_lables = []
    for (i, x) in enumerate(p1.series_list)
        if x[:label] ∈ known_lables
            p1.series_list[i][:label] = ""
        else
            push!(known_lables, x[:label])
        end
    end

    Plots.default(
        fontfamily = "Computer Modern",
        linewidth = 1,
        framestyle = :box,
        label = nothing,
        grid = false,
        legend = :topleft,
        margin=(Plots.Measures.Length(:mm, 2))
    )
    plot!(
        size = (600, 500),
        xlims = (LMIN  * 1.01, LMAX * 1.01),
    )

    savefig(string(@__DIR__) * "/plot.pdf") # Saves the plot from p as a .pdf vector graphic
end

main()