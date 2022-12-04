struct SIR{T<:AbstractFloat}
    β::T
    γ::T
end

function (model::SIR)(du, u, p, t)
    S, I, R = u
    (; β, γ) = model
        
    du[1] = -β * S * I
    du[2] = β * S * I - γ * I
    du[3] = γ * I

    return nothing
end

struct SIRV{T<:AbstractFloat}
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
