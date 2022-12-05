"""
    AbstractSIRVModel{T} <: AbstractCompartmentalModel{T}

Supertype of all SIRV models.
"""
abstract type AbstractSIRVModel{T} <: AbstractCompartmentalModel{T} end

@doc raw"""
    SIRV{T} <: AbstractSIRVModel{T}

SIRV model with contact rate ``\beta``, recovery rate ``\gamma``, vaccination
rate ``\nu``, and no demography:

```math
\begin{aligned}
\frac{dS}{dt} & = -\beta I S -\nu S\\[8pt]
\frac{dI}{dt} & = \beta I S - \gamma I \\[8pt]
\frac{dR}{dt} & = \gamma I \\[8pt]
\frac{dV}{dt} & = \nu S .
\end{aligned}
```

Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
struct SIRV{T} <: AbstractSIRVModel{T}
    β::T
    γ::T
    ν::T
end

@doc raw""" 
    SIRV(ν) -> SIRV

Return an `SIRV` model with vaccination rate ``\nu`` and default parameters
``\beta = 0.35`` and  ``\gamma = 0.035``.
"""
function SIRV(ν)
    return SIRV(0.35, 0.035, ν)
end

"""
    (::SIRV)(du, u, p, t) -> nothing

Update the derivates `du` given a model of type `SIRV`.
"""
function (model::SIRV)(du, u, p, t)
    S, I, R, V = u
    (; β, γ, ν) = model
        
    du[1] = -β * S * I - ν * S
    du[2] = β * S * I - γ * I
    du[3] = γ * I
    du[4] = ν * S

    return nothing
end

@doc raw"""
    SIRVSeasonalContact{T} <: AbstractSIRVModel{T}

SIRV model with seasonal contact rate ``\beta(t)``, recovery rate ``\gamma``, 
vaccination rate ``\nu``, and no demography.

The contact rate ``\beta(t)`` is given by 
``\beta(t) = \beta_0(1 + \beta_1 cos(2\pi t / 365))``, where ``\beta_0`` is the
mean contact rate and ``0 \leq \beta_1 \leq 1`` controls the degree of 
seasonality. 

The full equations of motion are:

```math
\begin{aligned}
\frac{dS}{dt} & = -\beta(t) I S -\nu S\\[8pt]
\frac{dI}{dt} & = \beta(t) I S - \gamma I \\[8pt]
\frac{dR}{dt} & = \gamma I \\[8pt]
\frac{dV}{dt} & = \nu S .
\end{aligned}
```

Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
struct SIRVSeasonalContact{T} <: AbstractSIRVModel{T}
    β₀::T
    β₁::T
    γ::T
    ν::T
end

@doc raw""" 
    SIRVSeasonalContact(β₁, ν) -> SIRVSeasonalContact

Return an `SIRVSeasonalContact` model with seasonality parameter ``\beta_1``, 
vaccination rate ``\nu``, and default parameters ``\beta_0 = 0.35`` and 
``\gamma = 0.035``.
"""
function SIRVSeasonalContact(β₁, ν)
    return SIRVSeasonalContact(0.35, β₁, 0.035, ν)
end

"""
    (::SIRVSeasonalContact)(du, u, p, t) -> nothing

Update the derivates `du` given a model of type `SIRVSeasonalContact`.
"""
function (model::SIRVSeasonalContact)(du, u, p, t)
    S, I, R, V = u
    (; β₀, β₁, γ, ν) = model

    β = β₀ * (1 + β₁ * cos(2π * t / 365))
        
    du[1] = -β * S * I - ν * S
    du[2] = β * S * I - γ * I
    du[3] = γ * I
    du[4] = ν * S

    return nothing
end

@doc raw"""
    SIRVDecayingImmunity{T} <: SIRVDecayingImmunity{T}

SIRV model with contact rate ``\beta``, recovery rate ``\gamma``, vaccination 
rate ``\nu``, infection immunity decay rate ``\alpha``, vaccination immunity 
decay rate ``\mu``, and no demography:

```math
\begin{aligned}
\frac{dS}{dt} & = -\beta I S -\nu S + \alpha R + \mu V\\[8pt]
\frac{dI}{dt} & = \beta I S - \gamma I \\[8pt]
\frac{dR}{dt} & = \gamma I - \alpha R\\[8pt]
\frac{dV}{dt} & = \nu S - \mu V.
\end{aligned}
```

Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
struct SIRVDecayingImmunity{T} <: AbstractSIRVModel{T}
    β::T
    γ::T
    ν::T
    α::T
    μ::T
end

@doc raw""" 
    SIRVDecayingImmunity(ν, α, μ) -> SIRVDecayingImmunity

Return an `SIRVDecayingImmunity` model with vaccination rate ``\nu``, infection 
immunity decay rate ``\alpha``, vaccination immunity decay rate ``\mu``, and 
default parameters ``\beta_0 = 0.35`` and ``\gamma = 0.035``.
"""
function SIRVDecayingImmunity(ν, α, μ)
    return SIRVDecayingImmunity(0.35, 0.035, ν, α, μ)
end

"""
    (::SIRVDecayingImmunity)(du, u, p, t) -> nothing

Update the derivates `du` given a model of type `SIRVDecayingImmunity`.
"""
function (model::SIRVDecayingImmunity)(du, u, p, t)
    S, I, R, V = u
    (; β, γ, ν, α, μ) = model

    du[1] = -β * S * I - ν * S + α * R + μ * V
    du[2] = β * S * I - γ * I
    du[3] = γ * I - α * R
    du[4] = ν * S - μ * V

    return nothing
end

@doc raw"""
    SIRVSeasonalContactDecayingImmunity{T} <: SIRVSeasonalContactDecayingImmunity{T}

SIRV model with seasonal contact rate ``\beta(t)``, recovery rate ``\gamma``, 
vaccination rate ``\nu``, infection immunity decay rate ``\alpha``, vaccination
immunity decay rate ``\mu``, and no demography.

The contact rate ``\beta(t)`` is given by 
``\beta(t) = \beta_0(1 + \beta_1 cos(2\pi t / 365))``, where ``\beta_0`` is the
mean contact rate and ``0 \leq \beta_1 \leq 1`` controls the degree of 
seasonality. 

The full equations of motion are:

```math
\begin{aligned}
\frac{dS}{dt} & = -\beta(t) I S -\nu S + \alpha R + \mu V\\[8pt]
\frac{dI}{dt} & = \beta(t) I S - \gamma I \\[8pt]
\frac{dR}{dt} & = \gamma I - \alpha R\\[8pt]
\frac{dV}{dt} & = \nu S - \alpha V.
\end{aligned}
```

Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
struct SIRVSeasonalContactDecayingImmunity{T} <: AbstractSIRVModel{T}
    β₀::T
    β₁::T
    γ::T
    ν::T
    α::T
    μ::T
end

@doc raw""" 
    SIRVSeasonalContactDecayingImmunity(β₁, ν, α, μ) -> SIRVSeasonalContactDecayingImmunity

Return an `SIRVSeasonalContactDecayingImmunity` model with seasonality parameter
``\beta_1``, vaccination rate ``\nu``, infection immunity decay rate ``\alpha``,
vaccination immunity decay rate ``\mu``, and default parameters 
``\beta_0 = 0.35`` and ``\gamma = 0.035``.
"""
function SIRVSeasonalContactDecayingImmunity(β₁, ν, α, μ)
    return SIRVSeasonalContactDecayingImmunity(0.35, β₁, 0.035, ν, α, μ)
end

"""
    (::SIRVSeasonalContactDecayingImmunity)(du, u, p, t) -> nothing

Update the derivates `du` given a model of type `SIRVSeasonalContactDecayingImmunity`.
"""
function (model::SIRVSeasonalContactDecayingImmunity)(du, u, p, t)
    S, I, R, V = u
    (; β₀, β₁, γ, ν, α, μ) = model

    β = β₀ * (1 + β₁ * cos(2π * t / 365))

    du[1] = -β * S * I - ν * S + α * R + μ * V
    du[2] = β * S * I - γ * I
    du[3] = γ * I - α * R
    du[4] = ν * S - μ * V

    return nothing
end
