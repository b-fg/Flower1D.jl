using Plots

# Plot solution points continously
function plotFlow(uδD :: Array{T, 2} where {T}, mesh :: Mesh; kwargs...)
    u = zeros(eltype(uδD), meshDOFNoGhost(mesh))
    x = zeros(eltype(uδD), meshDOFNoGhost(mesh))
    k = meshOrder(mesh) + 1 # number of solution points per element

    for i in 2:meshSize(mesh) - 1
        u[(i - 2) * k + 1:(i - 2) * k + k] = uδD[i, :]
        x[(i - 2) * k + 1:(i - 2) * k + k] = mesh.elements[i].innerPoints[:]
    end

    plot!(x, u; label = "", kwargs...)
end

# Plot flow at arbitrary coirdinates using LagrangeInterpolationAt
function plotFlowAt(xi :: Vector, uδD :: Array{T, 2} where {T}, mesh :: Mesh; kwargs...)
    u = zeros(eltype(uδD), length(xi))
    xj = mesh.standardElement.innerPoints
    elementNodes = mesh.elementsNodes
    wj = lagrangeBarycentricWeights(xj)

    for (i, x) ∈ enumerate(xi)
        elementIndex = reduce(vcat,transpose.([elementNodes[:, 1] .<= x] + [elementNodes[:, 2] .>= x]))
        n = findall(x->x == 2, elementIndex)[1][2] # convert to actual element index
        ξ = ξFromx(x, elementNodes[n, 1], elementNodes[n, 2])
        u[i] = lagrangeInterpolationAt(ξ, xj, uδD[n + 1, :], wj)
    end

    plot!(xi, u; label = "", kwargs...)
end