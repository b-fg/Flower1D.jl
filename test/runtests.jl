using Flower1D, Test

@testset "standard element" begin
    standardEl = StandardElement(1, 2, "legendre")
    @test 1==1
end


# @testset "lagrange inteprolation" begin
#     using Plots, Polynomials

#     f(x) = 1 / (1 + 16 * x^2)

#     # Define constants
#     TF, TI = Float64, Int64 # define precision
#     D = 1 # dimensions
#     E = 5 # number of elements per dimension
#     GP = 2^D # number of geometrical nodes per element (line,quad,hexa)
#     ELTYPE = 1 # element type, (1: line, 2: quad, 3: triangle)
#     EORDER = 2 # element order
#     LMIN, LMAX = -1.0, 1.0
#     h = (LMAX - LMIN) / E

#     # Create mesh
#     nodes = collect(LinRange(LMIN, LMAX, E + 1))
#     elnodes = Array{TF,3}(undef, E, GP, D)
#     for i ∈ 1:length(nodes)-1
#         elnodes[i, :, :] = [nodes[i] nodes[i+1]]
#     end
#     eltype = repeat([ELTYPE], E)
#     elorder = repeat([EORDER], E)
#     m = Flower1D.Mesh(elnodes, eltype, elorder, "legendre"; TF = TF, TI = TI)

#     # Plot target function
#     n_test_points_f = 200
#     x = collect(LinRange(LMIN, LMAX, n_test_points_f))
#     y_=f.(x)
#     plot(x,y_,label="f(x)", color="black")

#     # Discrete element solution
#     n_test_points_e = trunc(Int, n_test_points_f/E) # test points for each element

#     for (i,e) in enumerate(m.elements)
#         local xi = e.spoints # solution points
#         local yi = f.(xi) # solution values

#         local x = collect(LinRange(e.gpoints[1], e.gpoints[end], n_test_points_e))
#         local y = lagrange_1D(xi, yi, x)
#         plot!(x, y, label="Element $i", markershape=:circle)
#     end

#     plot!(xlims = (-1, 1), ylims = (0, 1.2))
#     savefig(string(@__DIR__)*"/plot.pdf") # Saves the plot from p as a .pdf vector graphic

#     @test 1==1
# end


# @testset "lagrange inteprolation" begin
#     @test 1==1
# end



# Polynomial Lagrange interpolation with basis functions
# li = lagrange_basis(xi)
# global L = Polynomials.Polynomial(0)  # Lagrange resulting interpolant L(x)
# for j in 1:length(xi)
#     global L += yi[j] * li[j]
# end
# y = [L(xi) for xi in x]
# plot!(x, y)


# xi = [0.,1.,2.,5.]
# yi = [2.,3.,12.,147.]

# xi = [-2.0/3.0,0.0,2.0/3.0]
# yi = [2.,3.,12.,147.]
# y = Interpolants.lagrange_1D(xi,yi,[3.])


# ToDO: Add more elements