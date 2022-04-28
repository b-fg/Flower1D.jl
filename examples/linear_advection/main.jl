using Revise, Flower1D
using BenchmarkTools, Plots, Printf, LaTeXStrings, LinearAlgebra

# Initial conditions
f1(x) = 1 / (1 + 16 * x ^ 2)
f2(x) = sin(pi * x)
f3(x) = exp(-20.0 * x ^ 2)
f4(x) = 0.25 + 0.5 * sin(pi * (2.0 * x - 1.0))

# Constants
T = Float64 # define precision

PDE = 1 # PDE flag: 1 ≡ Scalar, 2 ≡ Burgers, 3 ≡ Viscous Burgers, 4 ≡ Euler
BC = 0 # BC flag: 0 ≡ Periodic, 1 ≡ Outlet
LMIN, LMAX = -1, 1 # domain boundaries
N = 32 # number of elements
E = 1 # elements type
P = 2 # element order
C = 2 # Correction function parameter
IC = Function[f3]

U = 1.0 # characteristic velocity
L = LMAX - LMIN # characteristic length
CTU = L / U # convective time unit
TMAX = 5 * CTU # maximum simulation time [non-dim by CTU units]
TPRINT = 1 * CTU # print simulation time [non-dim by CTU units]
TPRINTEPS = 1e-3 # print simulation time reminder
CFL = 0.2

function main()
    # Init plot
    a = plot()

    # Create mesh
    vertexPoints = collect(LinRange(LMIN, LMAX, N + 1))
    elementsVertices = Matrix{T}(undef, N, 2)
    for i ∈ 1:length(vertexPoints) - 1
        elementsVertices[i, :] = [vertexPoints[i], vertexPoints[i + 1]]
    end
    mesh = Mesh(E, P, "legendre", elementsVertices, C; T = T)
    innerPoints = vcat([collect(transpose(mesh.elements[n].innerPoints)) for n ∈ 1:meshSize(mesh)]...)

    # Initialize flow
    flow = Flow(N, P, PDE, BC; T = T)
    init!(flow, IC, innerPoints)

    # Plot flow
    plotFlow(flow.uδD[:, :, 1], mesh, label = @sprintf("t = %.1f", flow.time[end] / CTU) |> latexstring)

    i = 1
    @time while flow.time[end] < TMAX
        # stepEuler!(flow, mesh, CFL)
        stepRK2SSP!(flow, mesh, CFL)
        # stepRK3SSPS3!(flow, mesh, CFL)
        i += 1
        @printf "%4i | δt = %.3e | t CTU = %.3e \n" i flow.δt[end] flow.time[end] / CTU

        # Display solution
        if flow.time[end] % TPRINT < flow.δt[end]
            plotFlow(flow.uδD[:, :, 1], mesh, label = @sprintf("t = %.1f", flow.time[end] / CTU) |> latexstring)
        end
    end

    # Format and save fig
    Plots.default(
        fontfamily = "Computer Modern",
        linewidth = 1,
        framestyle = :box,
        label = nothing,
        grid = false,
        margin=(Plots.Measures.Length(:mm, 2))
    )
    plot!(
        xlabel = L"x",
        ylabel = L"u",
        size = (600, 500),
        xlims = (LMIN, LMAX),
    )

    solutionLabel = @sprintf("N = %i, P = %i, C = %i", N, P, C) |> latexstring
    annotate!(LMIN + 0.05, 0.98, text(solutionLabel, :black, :left, 8))

    savefig(string(@__DIR__) * "/plot.pdf") # Saves the plot from p as a .pdf vector graphic
end

# Main
main()