"""
    Flow{T <: AbstractFloat}

Derived type storing flow related quantities.

# Fields
- `uδD :: Array{T, 3}`. Discontinuous conserved quantities. Stores values at each [element, integration point, equation index].
- `wδD :: Array{T, 3}`. Discontinuous primitive quantities. Stores values at each [element, integration point, equation index].
- `fδD :: Array{T, 3}`. Discontinuous advective flux quantities. Stores values at [element, integration point, equation index].
- `qδD :: Array{T, 3}`. Discontinuous viscous flux quantities. Stores values at [element, integration point, equation index].
- `rhs :: Array{T, 3}`. RHS quantities (du/dt = RHS). Stores values at each [element, integration point, equation index].
- `uδI :: Matrix{T}`. Interaction (aka common, numerical) flux of conserved quantities. Stores values at each [face, equation index].
- `fδI :: Matrix{T}`. Interaction advective flux. Stores values at each [face, equation index].
- `qδI :: Matrix{T}`. Interaction viscous flux. Stores values at each [face, equation index].
- `uK :: Array{T, 3}`. Left/Right conserved quantities at the faces. Stores values at each [face, left/right side, equation index].
- `fK :: Array{T, 3}`. Left/Right advective flux at the faces. Stores values at each [face, left/right side, equation index].
- `qK :: Array{T, 3}`. Left/Right viscous flux at the faces. Stores values at each [face, left/right side, equation index].
- `wBC :: Array{T, 3}`. Primitive quantities stored at the physical boundaries (acting as BCs).
    Stores values at each [left/right boundary element, integration point, equation index].
- `δt :: Vector{T}`. Stores the time step for each update of the equations.
- `time :: Vector{T}`. Stores the simulation time for each update of the equations.
- `PDE :: Integer`. Flag defining the type of equations to be solved:
    1 ≡ Linear advection equation, 2 ≡ Burgers equation, 3 ≡ Viscous Burgers equation, 4 ≡ Euler equations.
- `BC :: Integer`. Flag defining the type of boundary conditions:
    0 ≡ Periodic, 1 ≡ Outlet.
- `γ :: AbstractFloat`. Ratio of specific heats.
- `ν :: AbstractFloat`. Viscosity.

# Arguments
- `N :: Integer`. Number of elements (excluding ghost elements).
- `P :: Integer`. Order of the elements.
- `PDE :: Integer`. (As above).
- `BC :: Integer`. (As above).
- `γ :: AbstractFloat`. (As above).
- `ν :: AbstractFloat`. (As above).
- `T :: Type`. Float precision type.
"""
struct Flow{T <: AbstractFloat}
    uδD :: Array{T, 3}
    wδD :: Array{T, 3}
    fδD :: Array{T, 3}
    qδD :: Array{T, 3}
    rhs :: Array{T, 3}
    uδI :: Matrix{T}
    fδI :: Matrix{T}
    qδI :: Matrix{T}
    uK :: Array{T, 3}
    fK :: Array{T, 3}
    qK :: Array{T, 3}
    wBC :: Array{T, 3}
    δt :: Vector{T}
    time :: Vector{T}
    PDE :: Integer
    BC :: Integer
    γ :: AbstractFloat
    ν :: AbstractFloat

"""
    Flow(N :: Integer, P :: Integer, PDE :: Integer, BC :: Integer;
        γ :: AbstractFloat = 1.4, ν :: AbstractFloat = 0.0, T = Float64)

Initializes the `Flow` struct.
"""
    function Flow(N :: Integer, P :: Integer, PDE :: Integer, BC :: Integer;
        γ :: AbstractFloat = 1.4, ν :: AbstractFloat = 0.0, T = Float64)

        @assert PDE >= 1 || PDE <= 4 "PDE = $PDE, not implemented."
        @assert BC >= 0 || BC <= 1 "BC = $BC, not implemented."

        # Define number of equations according to PDE.
        NE = PDE < 4 ? 1 : 3

        uδD = zeros(T, (N + 2, P + 1, NE))
        wδD = zeros(T, (N + 2, P + 1, NE))
        qδD = zeros(T, (N + 2, P + 1, NE))
        wBC = zeros(T, (2, P + 1, NE))
        fδD = zeros(T, (N + 2, P + 1, NE))
        rhs = zeros(T, (N + 2, P + 1, NE))
        uδI = zeros(T, (N + 1, NE))
        fδI = zeros(T, (N + 1, NE))
        qδI = zeros(T, (N + 1, NE))
        uK = zeros(T, (N + 1, 2, NE))
        fK = zeros(T, (N + 1, 2, NE))
        qK = zeros(T, (N + 1, 2, NE))
        new{T}(uδD, wδD, fδD, qδD, rhs, uδI, fδI, qδI, uK, fK, qK, wBC, T[], T[0.0], PDE, BC, γ, ν)
    end
end

"""
    init!(flow :: Flow, f :: Vector{Function}, innerPoints :: Matrix{T} where T <: AbstractFloat)

Flow initialization function. Once a `Flow` struct is created, this function initializes the `flow.uδD` values
according to a `Vector`` of analytical functions `f` and the type of equations to be solved defined in `flow.PDE`.
The number of analytical functions in `f` has to match the size of the `Flow` number of equations.
"""
function init!(flow :: Flow, f :: Vector{Function}, innerPoints :: Matrix{T} where T <: AbstractFloat)
    @assert length(f) == size(flow.wδD, 3) "Number of initial analytical functions does not match the number of primitive variables."

    for (s, fi) in enumerate(f)
        @. flow.wδD[:, :, s] = fi(innerPoints)
    end

    flow.wBC[1, :, :] .= flow.wδD[2, :, :]
    flow.wBC[2, :, :] .= flow.wδD[end - 1, :, :]

    for m ∈ 1:size(flow.uδD, 2), n ∈ 2:size(flow.uδD, 1) - 1
        u = conservative(flow.wδD[n, m, :], flow.PDE)
        flow.uδD[n, m, :] .= u
    end
end

"""
    δt!(flow :: Flow,  δx :: Vector{T} where T<: AbstractFloat,  p :: Integer, CFL :: AbstractFloat = 0.1)

Computes a suitable time step according to the fastes wave speed within the solution and the CFL condition.
"""

# @fastmath function δt!(flow :: Flow,  δx :: Vector{T} where T<: AbstractFloat,  p :: Integer, CFL :: AbstractFloat = 0.1)
#     u = 2.0 * maximum(flow.uδD)
#     dx = maximum(@. δx / (p + 1))
#     δt = CFL * dx / u # This is very restrictive, instead should be looking elementwise: Δx[i] / u[i]
#     push!(flow.δt, δt)
# end

@fastmath function δt!(flow :: Flow,  δx :: Vector{T} where T<: AbstractFloat,  p :: Integer, CFL :: AbstractFloat = 0.1)
    if flow.PDE < 4 # Linear advection or Burger (Viscous) equation
        u = maximum(flow.uδD)
        δtAdv = minimum(CFL .* δx ./ (p + 1) ./ (2.0 .* u))
        if flow.ν ≉ 0.0
            δtVis = minimum(@. 0.5 * δx * δx / flow.ν)
            push!(flow.δt, 1.0 / (1.0 / δtAdv + 1.0 / δtVis))
        else
            push!(flow.δt, δtAdv)
        end
    else # Euler equations
        pressureDensityRatio = flow.wδD[2:end - 1, :, end] ./ flow.wδD[2:end - 1, :, 1]
        pRho = maximum(pressureDensityRatio)
        soundSpeed = sqrt(flow.γ * pRho) # local element sound speed
        u = maximum(flow.wδD[2:end - 1, :, 2]) # local element velocity
        S = max(abs(u + soundSpeed), abs(u), abs(u - soundSpeed))
        δtAdv = CFL * δx[1] / (2.0 * S) # δx[1] assumes all element have same size
        push!(flow.δt, δtAdv)
    end
end

"""
    solutionBoundaryConditions!(flow :: Flow)

Applies the boundary conditions to the discontinuous conserved quantities (`flow.uδD`) according to
the `flow.BC` flag.
"""
@fastmath function solutionBoundaryConditions!(flow :: Flow)
    if flow.BC == 0 # Periodic BC
        flow.uδD[1, :, :] = flow.uδD[end - 1, :, :]
        flow.uδD[end, :, :] = flow.uδD[2, :, :]
    else # Outlet BC
        @inbounds @simd for m ∈ 1:size(flow.uδD, 2)
            uBCL = conservative(flow.wBC[1, m, :], flow.PDE)
            uBCR = conservative(flow.wBC[2, m, :], flow.PDE)
            flow.uδD[1, m, :] .= uBCL
            flow.uδD[end, m, :] .= uBCR
        end
        flow.uK[1, 1, :] .= flow.uδD[1, end, :]
        flow.uK[end, 2, :] .= flow.uδD[end, 1, :]
    end
end

"""
    fluxesBoundaryConditions!(flow :: Flow)

Applies the boundary conditions to the discontinuous advection and viscous fluxes (`flow.fδD` and `flow.qδD`)
according to the `flow.BC` flag.
"""
@fastmath function fluxesBoundaryConditions!(flow :: Flow)
    if flow.BC == 1
        flow.fK[1, 1, :] .= flux(flow.wBC[1, end, :], flow.PDE)
        flow.fK[end, 2, :] .= flux(flow.wBC[2, 1, :], flow.PDE)
        if flow.ν ≉ 0.0
            flow.qK[1, 1, :] .= flow.wBC[1, end, :]
            flow.qK[end, 2, :] .= flow.wBC[2, 1, :]
        end
    end
end

"""
    reconstructAtFace!(flow :: Flow, basisFace :: Matrix{T} where T<: AbstractFloat)

Reconstruct the discontinuous conserved quantities (`flow.uδD`) at element faces, yielding `flow.uK`
"""
@fastmath function reconstructSolutionAtFace!(flow :: Flow, basisFace :: Matrix{T} where T<: AbstractFloat)
    @inbounds for s ∈ 1:size(flow.uδD, 3), c ∈ 2:size(flow.fδI, 1) - 1
        elementIndexL, elementIndexR = face2elements(c)
        flow.uK[c, 1, s] = @views sum(flow.uδD[elementIndexL, :, s] .* basisFace[2, :])
        flow.uK[c, 2, s] = @views sum(flow.uδD[elementIndexR, :, s] .* basisFace[1, :])
    end
    @inbounds @simd for s ∈ 1:size(flow.uδD, 3)
        flow.uK[1, 2, s] = @views sum(flow.uδD[2, :, s] .* basisFace[1, :])
        flow.uK[end, 1, s] = @views sum(flow.uδD[end - 1, :, s] .* basisFace[2, :])
    end
    #= Here we stil need to set left side of c=1 and right side of c=end.
    This is done just one time when setting the boundary conditions (if not periodic) =#
    if flow.BC == 0
        flow.uK[1, 1, :] .= flow.uK[end, 1, :]
        flow.uK[end, 2, :] .= flow.uK[1, 2, :]
    end
end

"""
    reconstructFluxesAtFace!(flow :: Flow, basisFace :: Matrix{T} where T<: AbstractFloat)

Reconstruct the discontinuous fluxes (`flow.fδD` and `flow.qδD``) at element faces, yielding `flow.fK` and `flow.qK`
"""
@fastmath function reconstructFluxesAtFace!(flow :: Flow, basisFace :: Matrix{T} where T<: AbstractFloat)
    # Advection flux
    @inbounds for s ∈ 1:size(flow.uδD, 3), c ∈ 2:size(flow.fδI, 1) - 1
        elementIndexL, elementIndexR = face2elements(c)
        flow.fK[c, 1, s] = @views sum(flow.fδD[elementIndexL, :, s] .* basisFace[2, :])
        flow.fK[c, 2, s] = @views sum(flow.fδD[elementIndexR, :, s] .* basisFace[1, :])
    end
    @inbounds @simd for s ∈ 1:size(flow.uδD, 3)
        flow.fK[1, 2, s] = @views sum(flow.fδD[2, :, s] .* basisFace[1, :])
        flow.fK[end, 1, s] = @views sum(flow.fδD[end - 1, :, s] .* basisFace[2, :])
    end
    if flow.BC == 0
        flow.fK[1, 1, :] .= flow.fK[end, 1, :]
        flow.fK[end, 2, :] .= flow.fK[1, 2, :]
    end
    # Diffusive flux
    if flow.ν ≉ 0.0
        @inbounds for s ∈ 1:size(flow.uδD, 3), c ∈ 2:size(flow.fδI, 1) - 1
            elementIndexL, elementIndexR = face2elements(c)
            flow.qK[c, 1, s] = @views sum(flow.qδD[elementIndexL, :, s] .* basisFace[2, :])
            flow.qK[c, 2, s] = @views sum(flow.qδD[elementIndexR, :, s] .* basisFace[1, :])
        end
        @inbounds @simd for s ∈ 1:size(flow.uδD, 3)
            flow.qK[1, 2, s] = @views sum(flow.qδD[2, :, s] .* basisFace[1, :])
            flow.qK[end, 1, s] = @views sum(flow.qδD[end - 1, :, s] .* basisFace[2, :])
        end
        if flow.BC == 0
            flow.qK[1, 1, :] .= flow.qK[end, 1, :]
            flow.qK[end, 2, :] .= flow.qK[1, 2, :]
        end
    end
end

"""
    viscousFlux!(flow :: Flow, basisFace :: Array{T, 2} where T <: AbstractFloat)

Computes the discontinuous viscous flux `flow.qδD`, where ``q = \\partial u / \\partial x``.
Note that the interaction flux `flow.uδI` is first computed, which is the solution ot the Riemann problem
at the face given the left and right solution `flow.uK`.
"""
@fastmath function viscousFlux!(flow :: Flow, dhL :: Vector{T}, dhR :: Vector{T}, Dij :: Matrix{T}, JijInv :: Vector{T}) where T <: AbstractFloat
    # Compute interaction solution
    @inbounds @simd for c ∈ 1:size(flow.uδI, 1)
        flow.uδI[c, :] .= centralFlux.(flow.uK[c, 1, :], flow.uK[c, 2, :])
    end
    # Compute qδD
    @inbounds for s ∈ 1:size(flow.uδD, 3), n ∈ 2:size(flow.uδD, 1) - 1
        faceIndexL, faceIndexR = element2faces(n)
        # Jumps of conservative quantities at faces
        jumpL = flow.uδI[faceIndexL, s] - flow.uK[faceIndexL, 2, s]
        jumpR = flow.uδI[faceIndexR, s] - flow.uK[faceIndexR, 1, s]
        # Compute RHS for reference coordinates
        flow.qδD[n, :, s] .= Dij * flow.uδD[n, :, s] + jumpL .* dhL .+ jumpR .* dhR
        # Scale dynamic viscosity, and with Jij to obtain in physical coordinates
        flow.qδD[n, :, s] .*= -JijInv[n - 1] * flow.ν
    end
end

"""
    riemannFluxes!(flow :: Flow)

Computes the interaction (aka common, numerical) advective and viscous fluxes, ie `flow.fI` and `flow.qI`,
given the left and right state of these fluxes.
Note that `fK`, where `K` stands for left `L` or right `R`, are computed as ``f(uK)``, and not using the Lagrange polynomials at x=[-1, 1].
Eg: fL = f(-1) = sum(fj * lj(-1)), as explained in Huynh original paper.
Also requires entropy fix to mitigate Gibbs oscillations (not implemented yet).
"""
@fastmath function riemannFluxes!(flow :: Flow)
    @inbounds @simd for c ∈ 1:size(flow.fδI, 1)
        uL = flow.uK[c, 1, :]
        uR = flow.uK[c, 2, :]
        wL = primitive(uL, flow.PDE, flow.γ)
        wR = primitive(uR, flow.PDE, flow.γ)
        fL = flux(wL, flow.PDE)
        fR = flux(wR, flow.PDE)
        qL = flow.qK[c, 1, :]
        qR = flow.qK[c, 2, :]

        if flow.PDE == 1 # Linear advection
            @. flow.fδI[c, :] = roeFlux(uL, uR, fL, fR, 1.0)
        elseif flow.PDE == 2 # Burgers equation
            @. flow.fδI[c, :] = roeFlux(uL, uR, fL, fR, uL[1])
        elseif flow.PDE == 3 # Burgers viscous equation
            @. flow.fδI[c, :] = roeFlux(uL, uR, fL, fR, uL[1])
            @. flow.qδI[c, :] = centralFlux(qL, qR)
        else # Euler equations
            # flow.fδI[c, :] .= rusanov(uL, uR, fL, fR, wL, wR, flow.γ)
            # flow.fδI[c, :] .= hll(uL, uR, fL, fR, wL, wR, flow.γ)
            flow.fδI[c, :] .= hllc(uL, uR, fL, fR, wL, wR, flow.γ)
        end
    end
end

"""
    rhs!(flow :: Flow, mesh :: Mesh)

Compute the RHS of the equations so that the solution can be advanced using a time integration methods such as Euler or Runge-Kutta.
The multiple steps are defined next:
1. Apply the boundary conditions to the discontinuous solution `flow.uδD`.
2. Reconstruct discontinuous solution `flow.uδD` at faces, obtaining `flow.uK`.
3. Compute the viscous flux `flow.qδD`.
4. Compute the advection flux `flow.fδD`.
5. Reconstruct advection and viscous fluxes at faces, obtaining `flow.fK` and `flow.qK`.
6. Compute common fluxes at element faces `flow.fδI` and `flow.qδI`.
7. Compute the equation RHS, du/dt = RHS = -1/Jij * df/dr, where `Jij` is the Jacobian of the mesh.
"""
@fastmath function rhs!(flow :: Flow, mesh :: Mesh)
    # Get correction function derivative and the Lagrange derivative matrix at element integration points.
    basisFace = mesh.standardElement.basisFace
    dhL = mesh.polynomialData.correctionDerivativeL
    dhR = mesh.polynomialData.correctionDerivativeR
    Dij = mesh.polynomialData.Dij
    JijInv = mesh.JijInv

    # 1. Apply BCs to solution: flow.uδD
    solutionBoundaryConditions!(flow)

    # 2. Extrapolate conserved quantities at faces: flow.uK
    reconstructSolutionAtFace!(flow, basisFace)

    # 3. Compute viscous flux at solution points: flow.qδD
    if flow.ν ≉ 0.0 viscousFlux!(flow, dhL, dhR, Dij, JijInv) end

    # 4. Compute advection flux at solution points: flow.fδD
    @inbounds for m ∈ 1:size(flow.uδD, 2), n ∈ 1:size(flow.uδD, 1)
        flow.wδD[n, m, :] .= primitive(flow.uδD[n, m, :], flow.PDE)
        flow.fδD[n, m, :] .= flux(flow.wδD[n, m, :], flow.PDE)
    end

    # 5. Extrapolate advection and diffusion fluxes at faces: flow.fK, flow.qK
    reconstructFluxesAtFace!(flow, basisFace)

    # 6. Compute common advection and diffusion fluxes at faces: flow.fδI, flow.qδI
    riemannFluxes!(flow)

    # 7. Compute the equation RHS
    @inbounds for s ∈ 1:size(flow.uδD, 3), n ∈ 2:size(flow.uδD, 1) - 1
        faceIndexL, faceIndexR = element2faces(n)
        # Advection flux
        jumpLAdv = flow.fδI[faceIndexL, s] - flow.fK[faceIndexL, 2, s]
        jumpRAdv = flow.fδI[faceIndexR, s] - flow.fK[faceIndexR, 1, s]
        # Compute RHS for reference coordinates
        flow.rhs[n, :, s] .= Dij * flow.fδD[n, :, s] + jumpLAdv .* dhL .+ jumpRAdv .* dhR
        # Diffusion flux
        if flow.ν ≉ 0.0
            jumpLDiff = flow.qδI[faceIndexL, s] - flow.qK[faceIndexL, 2, s]
            jumpRDiff = flow.qδI[faceIndexR, s] - flow.qK[faceIndexR, 1, s]
            # Add diffusion flux to advection flux
            flow.rhs[n, :, s] .+= Dij * flow.qδD[n, :, s] + jumpLDiff .* dhL .+ jumpRDiff .* dhR
        end
        # Scale with Jij to obtain in physical coordinates
        flow.rhs[n, :, s] .*= -JijInv[n - 1]
    end
end

"""
    stepRK2SSP!(flow :: Flow, mesh :: Mesh, CFL :: AbstractFloat = 0.1)

Update the numerical solution `flow.uδD` using the Runge-Kutta strong stability preserving (SSP) 2nd order scheme.
The time step is computed using the CFL condition `CFL :: AbstractFloat`.
"""
@fastmath function stepRK2SSP!(flow :: Flow, mesh :: Mesh, CFL :: AbstractFloat = 0.1)
    u0 = flow.uδD[2:end - 1, :, :]

    rhs!(flow, mesh)
    δt!(flow, meshΔx(mesh), meshOrder(mesh), CFL)

    u1 = @. flow.uδD[2:end - 1, :, :] + flow.δt[end] * flow.rhs[2:end - 1, :, :]

    @. flow.uδD[2:end - 1, :, :] = u1
    rhs!(flow, mesh)

    u2 = @. u1 + flow.δt[end] * flow.rhs[2:end - 1, :, :]

    @. flow.uδD[2:end - 1, :, :] = 0.5 * (u0 + u2)
    push!(flow.time, flow.time[end] + flow.δt[end])
end

"""
    stepRK3SSPS3!(flow :: Flow, mesh :: Mesh, CFL :: AbstractFloat = 0.1)

Update the numerical solution `flow.uδD` using the Runge-Kutta strong stability preserving (SSP) 3rd order (3 stages) scheme.
The time step is computed using the CFL condition `CFL :: AbstractFloat`.
"""
@fastmath function stepRK3SSPS3!(flow :: Flow, mesh :: Mesh, CFL :: AbstractFloat = 0.1)
    u0 = flow.uδD[2:end - 1, :, :]

    rhs!(flow, mesh)
    δt!(flow, meshΔx(mesh), meshOrder(mesh), CFL)

    u1 = @. flow.uδD[2:end - 1, :, :] + flow.δt[end] * flow.rhs[2:end - 1, :, :]

    @. flow.uδD[2:end - 1, :, :] = u1
    rhs!(flow, mesh)

    u2 = @. 0.75 * u0 + 0.25 * (u1 + flow.δt[end] * flow.rhs[2:end - 1, :, :])

    @. flow.uδD[2:end - 1, :, :] = u2
    rhs!(flow, mesh)

    @. flow.uδD[2:end - 1, :, :] = 1.0 / 3.0 * u0 + 2.0 / 3.0 * (u2  + flow.δt[end] * flow.rhs[2:end - 1, :, :])
    push!(flow.time, flow.time[end] + flow.δt[end])
end

"""
    stepEuler!(flow :: Flow, mesh :: Mesh, CFL :: AbstractFloat = 0.1)

Update the numerical solution `flow.uδD` using the Euler scheme.
The time step is computed using the CFL condition `CFL :: AbstractFloat`.
"""
@fastmath function stepEuler!(flow :: Flow, mesh :: Mesh, CFL :: AbstractFloat = 0.1)
    δt!(flow, meshΔx(mesh), meshOrder(mesh), CFL)
    rhs!(flow, mesh)
    @. flow.uδD[2:end - 1, :, :] = flow.uδD[2:end - 1, :, :] + flow.δt[end] * flow.rhs[2:end - 1, :, :]
    push!(flow.time, flow.time[end] + flow.δt[end])
end

"""
    flattenSolution(u :: Array{T, 2} where T <: AbstractFloat)

Returns a `Vector` type containing flow.uδD with all integration points of the mesh flattened in a single `Vector`.
"""
@fastmath function flattenSolution(u :: Array{T, 2} where T <: AbstractFloat)
    uFlat = zeros(eltype(u), (size(u, 1) - 2) * size(u, 2))
    k = size(u, 2) # number of solution points per element
    @inbounds  for n ∈ 2:size(u, 1) - 1
        uFlat[(n - 2) * k + 1:(n - 2) * k + k] = u[n, :]
    end
    return uFlat
end

"""
    function updatePrimitiveVars!(flow :: Flow)

Update the primitive variables array of the a `Flow`, ie `Flow.wδD`.
"""
@fastmath function updatePrimitiveVars!(flow :: Flow)
    @inbounds for m ∈ 1:size(flow.uδD, 2), n ∈ 2:size(flow.uδD, 1) - 1
        flow.wδD[n, m, :] .= primitive(flow.uδD[n, m, :], flow.PDE)
    end
end