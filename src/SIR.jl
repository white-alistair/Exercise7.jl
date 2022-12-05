"""
    AbstractSIRModel{T} <: AbstractCompartmentalModel{T}

Supertype of all SIR models.
"""
abstract type AbstractSIRModel{T} <: AbstractCompartmentalModel{T} end

@doc raw"""
    SIR{T} <: AbstractSIRModel{T}

Standard SIR model with contact rate ``\beta``, recovery rate ``\gamma``, and no
demography:

```math
\begin{aligned}
\frac{dS}{dt} & = -\beta I S \\[8pt]
\frac{dI}{dt} & = \beta I S - \gamma I \\[8pt]
\frac{dR}{dt} & = \gamma I .
\end{aligned}
```

Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
struct SIR{T} <: AbstractSIRModel{T}
    β::T
    γ::T
end

@doc raw"""
    SIR() -> SIR

Return an `SIR` model with default parameters ``\beta = 0.35`` and 
``\gamma = 0.035``.
"""
function SIR()
    return SIR(0.35, 0.035)
end

"""
    (::SIR)(du, u, p, t) -> nothing

Update the derivates `du` given a model of type `SIR`.
"""
function (model::SIR)(du, u, p, t)
    S, I, R = u
    (; β, γ) = model

    du[1] = -β * S * I
    du[2] = β * S * I - γ * I
    du[3] = γ * I

    return nothing
end

@doc raw"""
    SIRDecayingImmunity{T} <: AbstractSIRModel{T}

SIR model with contact rate ``\beta``, recovery rate ``\gamma``, infection 
immunity decay rate ``\alpha``, and no demography:

```math
\begin{aligned}
\frac{dS}{dt} & = -\beta I S + \alpha R \\[8pt]
\frac{dI}{dt} & = \beta I S - \gamma I \\[8pt]
\frac{dR}{dt} & = \gamma I - \alpha R .
\end{aligned}
```

Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
struct SIRDecayingImmunity{T} <: AbstractSIRModel{T}
    β::T
    γ::T
    α::T
end

@doc raw""" 
    SIRDecayingImmunity(α) -> SIRDecayingImmunity

Return an `SIRDecayingImmunity` model with immunity decay rate ``\alpha`` and
default parameters ``\beta = 0.35`` and ``\gamma = 0.035``.
"""
function SIRDecayingImmunity(α)
    return SIRDecayingImmunity(0.35, 0.035, α)
end

"""
    (::SIRDecayingImmunity)(du, u, p, t) -> nothing

Update the derivates `du` given a model of type `SIRDecayingImmunity`.
"""
function (model::SIRDecayingImmunity)(du, u, p, t)
    S, I, R = u
    (; β, γ, α) = model

    du[1] = -β * S * I + α * R
    du[2] = β * S * I - γ * I
    du[3] = γ * I - α * R

    return nothing
end
