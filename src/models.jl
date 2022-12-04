abstract type CompartmentalModel{T<:AbstractFloat} end

struct SIR{T} <: CompartmentalModel{T}
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

struct SIRV{T} <: CompartmentalModel{T}
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
