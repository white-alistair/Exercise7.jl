module Exercise7

export SIR,
    SIRDecayingImmunity,
    SIRV,
    SIRVSeasonalContact,
    SIRVDecayingImmunity,
    SIRVSeasonalContactDecayingImmunity,
    plot,
    plot_phase_diagram,
    plot_trajectory_given_vax_rate,
    plot_total_infections_by_vax_rate

using DynamicalSystems, DifferentialEquations, CairoMakie, Printf

"""
    AbstractCompartmentalModel{T<:AbstractFloat}

Supertype of all compartmental models defined in this package.
"""
abstract type AbstractCompartmentalModel{T<:AbstractFloat} end

include("SIR.jl")
include("SIRV.jl")
include("plots.jl")

end
