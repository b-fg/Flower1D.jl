using Revise, Flower1D
using BenchmarkTools, Plots, Printf, LaTeXStrings, LinearAlgebra

gr()

# Constants
T = Float64 # define precision

PDE = 4 # PDE flag: 1 ≡ Scalar, 2 ≡ Burgers, 3 ≡ Viscous Burgers, 4 ≡ Euler
BC = 1 # BC flag: 0 ≡ Periodic, 1 ≡ Outlet
LMIN, LMAX = 0.0, 1.0 # domain boundaries
N = 64 # number of elements (use an even number, so that the initial discontinuity is located at a face)
E = 1 # elements type
P = 3 # element order
C = 1 # Correction function parameter

U = 1.0 # characteristic velocity
L = LMAX - LMIN # characteristic length
TMAX = 0.2 # maximum simulation time
CFL = 0.25

# Sod shock tube initial condition
density(x) = x <= 0.5 ? 1.0 : 0.125
velocity(x) = 0.0
pressure(x) = x <= 0.5 ? 1.0 : 0.1
IC = Function[density, velocity, pressure]

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

    # Time loop
    i = 0
    @time while flow.time[end] < TMAX
        stepRK2SSP!(flow, mesh, T(CFL))
        i += 1
        @printf "%4i | δt = %.3e | t = %.3e \n" i flow.δt[end] flow.time[end]
    end

    # Update primitive variables
    updatePrimitiveVars!(flow)

    # Compute, plot, calculate error using the analytical solution for density
    problem = ShockTubeProblem(
        geometry = (0.0, 1.0, 0.5),
        left_state = (ρ = 1.0, u = 0.0, p = 1.0),
        right_state = (ρ = 0.125, u = 0.0, p = 0.1),
        t = TMAX,
        γ = 1.4)
    xi = integrationPointsLocation(mesh)
    positions, regions, values = solveSodShockTube(problem, xi)
    ρErr = abs.(values[:ρ] .- flattenSolution(flow.wδD[:, :, 1]))
    error = LinearAlgebra.norm(ρErr)
    @printf "L2 error = %.3e \n" error

    plot!(xi, values[:ρ]; label = "Ref @ " * (@sprintf("t = %.2f", flow.time[end]) |> latexstring), color = :black)
    vline!([positions["Shock"], positions["Foot of rarefaction"], positions["Head of rarefaction"],
        positions["Contact Discontinuity"]], color = :grey, linewidth = 1, linestyle = :dash)

    # Plot numerical solution for density and error
    solutionLabel = @sprintf("N = %i, P = %i, C = %i", N, P, C) |> latexstring
    plotFlow(flow.uδD[:, :, 1], mesh; color = :orange, linewidth = 2, label = solutionLabel)
    # plotFlowAt(collect(range(LMIN, LMAX, 200)), flow.uδD[:, :, 1], mesh; markershape = :circle, markersize = 3, label = solutionLabel)
    plot!(xi, ρErr; label = "L2 error", color = :red, linewidth = 1)

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
        xlabel = L"x",
        ylabel = L"\rho",
        size = (600, 500),
        xlims = (LMIN, LMAX),
        ylims = (0.0, 1.05),
    )
    savefig(string(@__DIR__) * "/plot.pdf")
end

# Main
main()