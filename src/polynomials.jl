"""
    lagrangeBasisAt(x :: AbstractFloat, xi :: Vector{T} where T <: AbstractFloat)

Lagrange basis polynomials values evaluated at location `AbstractFloat :: x` given `Vector :: xi` integration points.
https://en.wikipedia.org/wiki/Lagrange_polynomial
"""
function lagrangeBasisAt(x :: AbstractFloat, xi :: Vector{T} where T <: AbstractFloat)
    n = length(xi) # number of integration points
    lj = ones(eltype(xi),n)
    for i ∈ 1:n, j ∈ 1:n
        if i != j
            lj[i] *= (x - xi[j]) / (xi[i] - xi[j])
        end
    end
    return lj
end

"""
    dLagrangeBasisAt(x :: AbstractFloat, xj :: Vector{T} where T <: AbstractFloat)

Derivative of the Lagrange basis polynomials evaluated at location `x :: AbstractFloat` given `xj :: Vector` integration points.
"""
function dLagrangeBasisAt(x :: AbstractFloat, xj :: Vector{T} where T <: AbstractFloat)
    k = length(xj) # number of integration points
    dlj = zeros(eltype(xj), k)
    lj = lagrangeBasisAt(x, xj)
    for j ∈ 1:k
        sum = 0.0
        for i ∈ 1:k
            if (i != j) sum += 1.0 / (x - xj[i]) end
        end
        dlj[j] = lj[j] * sum
    end
    return dlj
end

# Implementation of `dLagrangeBasisAt` for multiple points `x :: Vector{T}` where to evaluate the derivative.
function dLagrangeBasisAt(x :: Vector{T}, xj :: Vector{T}) where T <: AbstractFloat
    return mapreduce(permutedims, vcat, dLagrangeBasisAt.(Ref(xj), x))
end

"""
    lagrangeBarycentricWeights(xj :: Vector{T} where T <: AbstractFloat)

Returns de barycentric weights required by the Lagrange barycentric interpolation given `xj :: Vector{T}` node locations.
This function implements the *Algorithm 30* in _Implementing Spectral Methods for Partial Differential Equation_ by D. A. Kopriva.
"""
function lagrangeBarycentricWeights(xj :: Vector{T} where T <: AbstractFloat)
    N = length(xj) # number of nodes
    wj = ones(N)
    for j ∈ 2:N, k ∈ 1:j - 1
        wj[k] *=  xj[k] - xj[j]
        wj[j] *=  xj[j] - xj[k]
    end
    return 1.0 ./ wj
end

"""
    lagrangeInterpolationAt(x :: AbstractFloat, xj :: Vector{T}, fj :: Vector{T}, wj :: Vector{T}) where T <: AbstractFloat

Returns the Lagrange interpolant function value at location `x :: AbstractFloat` given the fitting values `fj :: Vector{T}`,
    their locations `xj :: Vector{T}` and the barycentric weights `wj :: Vector{T}` computed with `lagrangeBarycentricWeights`.
This function implements the *Algorithm 31* in _Implementing Spectral Methods for Partial Differential Equation_ by D. A. Kopriva.
"""
function lagrangeInterpolationAt(x :: AbstractFloat, xj :: Vector{T}, fj :: Vector{T}, wj :: Vector{T}) where T <: AbstractFloat
    N = length(xj) # number of nodes
    numerator = 0.0
    denominator = 0.0
    for j ∈ 1:N
        if x ≈ xj[j] return fj[j] end
        t = wj[j] / (x - xj[j])
        numerator += t * fj[j]
        denominator += t
    end
    return numerator / denominator
end

# Implementation of `lagrangeInterpolationAt` for multiple points `x :: Vector{T}` where to evaluate the interpolant.
function lagrangeInterpolationAt(x :: Vector{T}, xj :: Vector{T}, fj :: Vector{T}, wj :: Vector{T}) where T <: AbstractFloat
    return lagrangeInterpolationAt.(x, Ref(xj), Ref(fj), Ref(wj))
end

"""
    lagrangeDerivativeMatrix(xj :: Vector{T} where T <: AbstractFloat)

Returns the Lagrange basis derivative matrix given a set of node locations `xj :: Vector{T}`.
It can be used to efficiently compute the Lagrange interpolant derivative at a node location ``j \\in [1, N]``
with a matrix-vector multiplication, using the function node weights:
``f'(x_j) = \\sum_{j=0}^{N} D_{ij} f_j, i=0, 1, ..., N``.

This function implements the *Algorithm 37* in _Implementing Spectral Methods for Partial Differential Equation_ by D. A. Kopriva.
"""
function lagrangeDerivativeMatrix(xj :: Vector{T} where T <: AbstractFloat)
    N = length(xj)
    Dij = zeros(N, N)
    wj = lagrangeBarycentricWeights(xj)
    for i ∈ 1:N
        Dij[i, i] = 0.0
        for j ∈ 1:N
            if i != j
                Dij[i, j] = wj[j] / (wj[i] * (xj[i] - xj[j]))
                Dij[i, i] -=  Dij[i, j]
            end
        end
    end
    return Dij
end

"""
    pLegendreAt(p :: Int, x :: Real)

Returns the values of the `p :: Int` Legendre polynomial evaluated at `x :: AbstractFloat`.
"""
function pLegendreAt(p :: Integer, x :: AbstractFloat)
    Lp = 0.0
    for k ∈ 0:p
        Lp += binomial(p, k) * binomial(p + k, k) * (0.5 * (x - 1)) ^ k
    end
    return Lp
end

"""
    dLegendreAt(p :: Integer, x :: AbstractFloat)

Returns the values of the 1st derivative of the `p :: Int` Legendre polynomial evaluated at `x <: AbstractFloat`.
"""
function dLegendreAt(p :: Integer, x :: AbstractFloat)
    dLp = 0.0
    for k ∈ 0:p
        dLp += binomial(p, k) * binomial(p + k, k) * 0.5 * k * (0.5 * (x - 1)) ^ (k - 1)
    end
    return dLp
end

"""
    vcjh(p :: Integer, x :: AbstractFloat; C :: Union{Integer, AbstractFloat} = 1,
        derivativeChoice :: Bool = false)

VCJH correction functions required to compute correction flux fδC.
The function can return the correction functions `hL, hR` or their derivative `dhL, dhR` depending
on the `derivativeChoice :: Bool` flag. The correction scheme is governed by the Union type `C`.
The following schemes are recovered if:
    - C = 0 :: Integer ➡ DG nodal
    - C = 1 :: Integer ➡ Spectral Difference (SD)
    - C = 2 :: Integer ➡ Huynh
    - C <: AbstractFloat ➡ Use this exact value.

Other arguments: `p :: Integer` order of the solution polynomial within each element.
`x :: AbstractFloat` is the location where the Legendre polynomials are evaluated (inner points of elements).
"""
function vcjh(p :: Integer, x :: AbstractFloat; C :: Union{Integer, AbstractFloat} = 1,
    derivativeChoice :: Bool = false)
    ap = factorial(2 * p) / ((2 ^ p) * (factorial(p) ^ 2))

    # Select the scheme according to c
    if typeof(C) <: AbstractFloat
        c = C
    elseif typeof(C) <: Integer && C == 0
        c = 0.0
    elseif typeof(C) <: Integer && C == 1
        c = 2 * p / ((2 * p + 1) * (p + 1) * ((ap * factorial(p)) ^ 2))
    elseif typeof(C) <: Integer && C == 2
        c = 2 * (p + 1) / ((2 * p + 1) * p * ((ap * factorial(p)) ^ 2))
    else
        error("VCJH scheme not implemented.")
    end

    ηp = 0.5 * c * (2 * p + 1) * ((ap * factorial(p)) ^ 2)

    # Select Legendre polynomials or its 1stO derivatives
    legendreAt = derivativeChoice ? dLegendreAt : pLegendreAt

    # Compute functions
    Lp = legendreAt(p, x)
    Lpp = legendreAt(p + 1, x)
    Lpm = legendreAt(p - 1, x)

    hL = 0.5 * ((-1.0) ^ p) * (Lp - (ηp * Lpm + Lpp) / (1 + ηp))
    hR = 0.5 * (Lp + (ηp * Lpm + Lpp) / (1 + ηp))
    return [hL, hR]
end

# Implementation of `vcjh` for multiple points `x :: Vector{T}` where to evaluate the correction function.
function vcjh(p :: Integer, x :: Vector{T} where T <: AbstractFloat; Args...)
    return mapreduce(permutedims, vcat, vcjh.(p, x; Args...))
end

"""
Converts element geometrical coordinates x to reference coordinates ξ in the [-1,1] interval
given the elements geometrical nodes (xi, xe).
"""
ξFromx(x, xi, xe) = 2.0 * (x - xi) / (xe - xi) - 1.0

"""
Converts reference coordinates ξ defined in interval [-1,1] to physical coordinates x
given the elements geometrical nodes (xi, xe).
"""
xFromξ(ξ, xi, xe) = 0.5 * (1.0 - ξ) * xi + 0.5 * (1.0 + ξ) * xe