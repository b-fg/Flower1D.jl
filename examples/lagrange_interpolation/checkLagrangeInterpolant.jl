using Revise, Flower1D

# Polynomial Lagrange interpolation at a point `x`
function checkLagrangeInterpolant()
    xi = [0., 1., 2., 5.]
    yi = [2., 3., 12., 147.]

    wi = lagrangeBarycentricWeights(xi)
    x = 3.0
    return lagrangeInterpolationAt(x, xi, yi, wi)
end

checkLagrangeInterpolant()