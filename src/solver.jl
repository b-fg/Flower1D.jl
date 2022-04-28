"""
    flux(w :: Vector{T} where T <: AbstractFloat, PDE :: Integer, γ :: AbstractFloat = 1.4)

Defines the flux according to the primitive variables `w` and the `PDE` to be solved.
For the Euler equations, the primitive variables are defined as:
- w[1]: Density.
- w[2]: Velocity X.
- w[3]: Pressure.
"""
@fastmath function flux(w :: Vector{T} where T <: AbstractFloat, PDE :: Integer, γ :: AbstractFloat = 1.4)
    # Linear advection
    if PDE == 1
        return w
    # Burgers equation
    elseif PDE < 4
        return 0.5 .* w .* w
    # Euler equations
    else
        f = similar(w)
        f[1] = w[1] * w[2]
        f[2] = w[1] * w[2] * w[2] + w[end]
        E = w[3] / (γ - 1.0) + 0.5 * w[1] * w[2] * w[2]
        f[3] = (E + w[end]) * w[2]
        return f
    end
end

"""
    roeFlux(uL :: AbstractFloat, uR :: AbstractFloat, fL :: AbstractFloat, fR :: AbstractFloat, au :: AbstractFloat)

Roe solver.
"""
@fastmath function roeFlux(uL :: AbstractFloat, uR :: AbstractFloat, fL :: AbstractFloat, fR :: AbstractFloat, au :: AbstractFloat)
    a = !(uL ≈ uR) ? (fR - fL) / (uR - uL + eps(typeof(uL))) : au
    # flux splitting
    return 0.5 * (fL + fR) - 0.5 * abs(a) * (uR - uL)
end

"""
    centralFlux(fL :: AbstractFloat, fR :: AbstractFloat)

Riemann solver using the central flux approximation (only adequate for viscous fluxes).
"""
@fastmath function centralFlux(fL :: AbstractFloat, fR :: AbstractFloat)
    return 0.5 * (fL + fR)
end

"""
    rusanov(uL :: Vector{T}, uR :: Vector{T}, fL :: Vector{T}, fR :: Vector{T}, wL :: Vector{T}, wR :: Vector{T},
        γ :: AbstractFloat = 1.4) where {T <: AbstractFloat}

Rusanov solver for the Euler equations.
"""
@fastmath function rusanov(uL :: Vector{T}, uR :: Vector{T}, fL :: Vector{T}, fR :: Vector{T}, wL :: Vector{T}, wR :: Vector{T},
    γ :: AbstractFloat = 1.4) where {T <: AbstractFloat}
    # Estimate wave speed S
    aL = sqrt(γ * wL[end] / wL[1])
    aR = sqrt(γ * wR[end] / wR[1])
    S = max(abs(wL[2]) + aL, abs(wR[2]) + aR)
    return @. 0.5 * (fL + fR) - 0.5 * S * (uR - uL)
end

"""
    function hll(uL :: Vector{T}, uR :: Vector{T}, fL :: Vector{T}, fR :: Vector{T}, wL :: Vector{T}, wR :: Vector{T},
        γ :: AbstractFloat = 1.4) where {T <: AbstractFloat}

HLL solver for the Euler equations.
"""
@fastmath function hll(uL :: Vector{T}, uR :: Vector{T}, fL :: Vector{T}, fR :: Vector{T}, wL :: Vector{T}, wR :: Vector{T},
    γ :: AbstractFloat = 1.4) where {T <: AbstractFloat}
    # Estimates wave speeds SL & SR using Roe averages
    uRoe = (sqrt(wL[1]) * wL[2] + sqrt(wR[1]) * wR[2]) / (sqrt(wL[1]) + sqrt(wR[1]))
    HL = (uL[end] + wL[end]) / wL[1] # Enthalpy left
    HR = (uR[end] + wR[end]) / wR[1] # Enthalpy right
    HRoe = (sqrt(wL[1]) * HL + sqrt(wR[1]) * HR) / (sqrt(wL[1]) + sqrt(wR[1])) # Averaged enthalpy
    aRoe = sqrt((γ - 1.0) * (HRoe - 0.5 * uRoe * uRoe)) # Avreaged sound speed
    SL = uRoe - aRoe # Left wave speed
    SR = uRoe + aRoe # Right wave speed

    # Alternative way to compute the wave speeds.
    # aL = sqrt(γ * wL[end] / wL[1])
    # aR = sqrt(γ * wR[end] / wR[1])
    # SL = min(wL[2] - aL, wR[2] - aR)
    # SR = max(wL[2] + aL, wR[2] + aR)

    # Compute the HLL flux
    if 0.0 <= SL
        return fL
    elseif SL <= 0.0 <= SR
        @. return (SR * fL - SL * fR + SL * SR * (uR - uL)) / (SR - SL)
    else
        return fR
    end
end

"""
    function hllc(uL :: Vector{T}, uR :: Vector{T}, fL :: Vector{T}, fR :: Vector{T}, wL :: Vector{T}, wR :: Vector{T},
        γ :: AbstractFloat = 1.4) where {T <: AbstractFloat}

HLLC solver for the Euler equations.
"""
@fastmath function hllc(uL :: Vector{T}, uR :: Vector{T}, fL :: Vector{T}, fR :: Vector{T}, wL :: Vector{T}, wR :: Vector{T},
    γ :: AbstractFloat = 1.4) where {T <: AbstractFloat}
    # Estimates wave speeds SL & SR using Roe averages
    uRoe = (sqrt(wL[1]) * wL[2] + sqrt(wR[1]) * wR[2]) / (sqrt(wL[1]) + sqrt(wR[1]))
    HL = (uL[end] + wL[end]) / wL[1] # Enthalpy left
    HR = (uR[end] + wR[end]) / wR[1] # Enthalpy right
    HRoe = (sqrt(wL[1]) * HL + sqrt(wR[1]) * HR) / (sqrt(wL[1]) + sqrt(wR[1])) # Averaged enthalpy
    aRoe = sqrt((γ - 1.0) * (HRoe - 0.5 * uRoe * uRoe)) # Avreaged sound speed
    SL = uRoe - aRoe # Left wave speed
    SR = uRoe + aRoe # Right wave speed
    Sstar = (wR[end] - wL[end] + uL[2] * (SL - wL[2]) - uR[2] * (SR - wR[2])) / (wL[1] * (SL - wL[2]) - wR[1] * (SR - wR[2]))

    # Estimate intermediate states *
    uLstar = similar(uL)
    uLstar[1] = 1.0
    uLstar[2] = Sstar
    uLstar[end] = uL[end] / wL[1] + (Sstar - wL[2]) * (Sstar + wL[end] / (wL[1] * (SL - wL[2])))
    @. uLstar *= wL[1] * ((SL - wL[2]) / (SL - Sstar))
    uRstar = similar(uR)
    uRstar[1] = 1.0
    uRstar[2] = Sstar
    uRstar[end] = uR[end] / wR[1] + (Sstar - wR[2]) * (Sstar + wR[end] / (wR[1] * (SR - wR[2])))
    @. uRstar *= wR[1] * ((SR - wR[2]) / (SR - Sstar))

    fLstar = similar(fL)
    @. fLstar = fL + SL * (uLstar - uL)
    fRstar = similar(fR)
    @. fRstar = fR + SR * (uRstar - uR)

    # Compute the HLLC flux
    if 0.0 <= SL
        return fL
    elseif SL <= 0.0 <= Sstar
        return fLstar
    elseif Sstar <= 0.0 <= SR
        return fRstar
    else
        return fR
    end
end

"""
    primitive(u :: Vector{T} where T <: AbstractFloat, PDE :: Integer, γ :: AbstractFloat = 1.4)

Computes the primitive variables `w` given the conservative variables `u` according to the `PDE`.
"""
@fastmath function primitive(u :: Vector{T} where T <: AbstractFloat, PDE :: Integer, γ :: AbstractFloat = 1.4)
    # Linear or Burgers equation
    if PDE < 4
        return u
    # Euler equations
    else
        w = similar(u)
        # Density
        w[1] = u[1]
        # Velocities
        for i ∈ 2:length(u) - 1
            w[i] = u[i] / u[1]
        end
        # Pressure
        velocitySquared = sum(x -> x^2, @view w[2:end - 1])
        w[end] = (γ - 1.0) * (u[end] - 0.5 * u[1] * velocitySquared)
        return w
    end
end

"""
    conservative(w :: Vector{T} where T <: AbstractFloat, PDE :: Integer, γ :: AbstractFloat = 1.4)

Computes the conservative variables `u` given the primitive variables `w` according to the `PDE`.
"""
@fastmath function conservative(w :: Vector{T} where T <: AbstractFloat, PDE :: Integer, γ :: AbstractFloat = 1.4)
    # Linear or Burgers equation
    if PDE < 4
        return w
    # Euler equations
    else
        u = similar(w)
        # Mass
        u[1] = w[1]
        # Momentum
        for i ∈ 2:length(u) - 1
            u[i] = w[1] * w[i]
        end
        # Energy
        velocitySquared = sum(x -> x^2, @view w[2:end - 1])
        u[end] = w[end] / (γ - 1.0) + 0.5 * w[1] * velocitySquared
        return u
    end
end