using Revise, Flower1D
using BenchmarkTools, Plots, Printf, LaTeXStrings, LinearAlgebra

# Initial conditions
f1(x) = 1 / (1 + 16 * x ^ 2)
f2(x) = sin(pi * x)
f3(x) = exp(-20.0 * x ^ 2)
f4(x) = 0.25 + 0.5 * sin(pi * (2.0 * x - 1.0))

# Constants
T = Float64 # define precision

PDE = 3 # PDE flag: 1 ≡ Scalar, 2 ≡ Burgers, 3 ≡ Viscous Burgers, 4 ≡ Euler
ν = 0.5 # Viscosity
BC = 0 # BC flag: 0 ≡ Periodic, 1 ≡ Outlet
LMIN, LMAX = 0.0, 2.0 * pi # domain boundaries
N = 64 # number of elements
E = 1 # elements type
P = 2 # element order
C = 2 # Correction function parameter

TMAX = 0.4 # maximum simulation time
TPRINT = 0.1 # print simulation time
CFL = 0.2

# Analytical function
function f(ν, t, x)
    a = x - 4.0 * t
    b = x - 4.0 * t  - 2.0 * pi
    c = 4.0 * ν * (t + 1.0)
    ϕ = exp(-a * a / c) + exp(-b * b / c)
    dϕ = -2.0 * a / c * exp(-a * a / c) - 2.0 * b / c  * exp(-b * b / c)
    return 4.0 - 2.0 * ν * dϕ / ϕ
end

#Initial condition
curry2(f, x1, x2) = (args...) -> f(x1, x2, args...)
f0 = curry2(f, ν, 0.0)
IC = Function[f0]

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
    flow = Flow(N, P, PDE, BC; ν = ν, T = T)
    init!(flow, IC, innerPoints)

    # Plot flow
    plotFlow(flow.uδD[:, :, 1], mesh, label = @sprintf("t = %.1f", flow.time[end]) |> latexstring)

    i = 1
    @time while flow.time[end] < TMAX
        # stepEuler!(flow, mesh, CFL)
        stepRK2SSP!(flow, mesh, CFL)
        # stepRK3SSPS3!(flow, mesh, CFL)
        i += 1
        @printf "%4i | δt = %.3e | t = %.3e \n" i flow.δt[end] flow.time[end]

        # Display solution
        if flow.time[end] % TPRINT < flow.δt[end]
            plotFlow(flow.uδD[:, :, 1], mesh, label = @sprintf("t = %.1f", flow.time[end]) |> latexstring)
        end
    end

    # Compute, plot, calculate error using the analytical solution
    xi = integrationPointsLocation(mesh)
    uErr = abs.(f.(ν, flow.time[end], xi) .- flattenSolution(flow.uδD[:, :, 1]))
    plot!(xi, f.(ν, 0.0, xi); label = "Ref @ " * (@sprintf("t = %.2f", 0.0) |> latexstring), color = :black, linewidth = 1, linestyle = :dash)
    plot!(xi, f.(ν, flow.time[end], xi); label = "Ref @ " * (@sprintf("t = %.2f", flow.time[end]) |> latexstring),
        color = :black, linewidth = 1, linestyle = :dashdot)
    error = LinearAlgebra.norm(uErr)
    @printf "L2 error = %.3e \n" error

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
        xlabel = L"x",
        ylabel = L"u",
        size = (600, 500),
        xlims = (LMIN, LMAX),
        xticks = ([0:pi / 2:2 * pi;], pifrac.([0:pi / 2:2 * pi;])),
        legend = :bottomleft
    )

    solutionLabel = @sprintf("N = %i, P = %i, C = %i", N, P, C) |> latexstring
    annotate!(LMIN + 0.2, maximum(f.(ν, 0.0, xi)), text(solutionLabel, :black, :left, 8))
    annotate!(LMIN + 0.2, maximum(f.(ν, 0.0, xi)) - 0.15, text(@sprintf("\\nu = %.2f", flow.ν) |> latexstring, :black, :left, 8))

    savefig(string(@__DIR__) * "/plot.pdf") # Saves the plot from p as a .pdf vector graphic
end

# Main
main()