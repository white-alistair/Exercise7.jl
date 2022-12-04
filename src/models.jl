abstract type AbstractCompartmentalModel{T<:AbstractFloat} end

# SIR Models
abstract type AbstractSIRModel{T} <: AbstractCompartmentalModel{T} end

struct SIR{T} <: AbstractSIRModel{T}
    β::T
    γ::T
end

function SIR()
    return SIR(0.35, 0.035)
end

function (model::SIR)(du, u, p, t)
    S, I, R = u
    (; β, γ) = model
        
    du[1] = -β * S * I
    du[2] = β * S * I - γ * I
    du[3] = γ * I

    return nothing
end

# SIRV Models
abstract type AbstractSIRVModel{T} <: AbstractCompartmentalModel{T} end

struct SIRV{T} <: AbstractSIRVModel{T}
    β::T
    γ::T
    ν::T
end

function (model::SIRV)(du, u, p, t)
    S, I, R, V = u
    (; β, γ, ν) = model
        
    du[1] = -β * S * I - ν * S
    du[2] = β * S * I - γ * I
    du[3] = γ * I
    du[4] = ν * S

    return nothing
end

struct SIRVSeasonal{T} <: AbstractSIRVModel{T}
    β₀::T
    β₁::T
    γ::T
    ν::T
end

function (model::SIRVSeasonal)(du, u, p, t)
    S, I, R, V = u
    (; β₀, β₁, γ, ν) = model

    β = β₀ * (1 + β₁ * cos(2π * t / 365))
        
    du[1] = -β * S * I - ν * S
    du[2] = β * S * I - γ * I
    du[3] = γ * I
    du[4] = ν * S

    return nothing
end

struct SIRVDecayingImmunity{T} <: AbstractSIRVModel{T}
    β::T
    γ::T
    ν::T
    α::T
    μ::T
end

function (model::SIRVDecayingImmunity)(du, u, p, t)
    S, I, R, V = u
    (; β, γ, ν, α, μ) = model
        
    du[1] = -β * S * I - ν * S + α * R + μ * V
    du[2] = β * S * I - γ * I
    du[3] = γ * I - α * R
    du[4] = ν * S - μ * V

    return nothing
end
