module Exercise7 
    export SIR, SIRV, SIRVSeasonal, SIRVDecayingImmunity, plot, plot_trajectory_given_vax_rate, plot_total_infections_by_vax_rate

    using DynamicalSystems, CairoMakie, Printf

    include("models.jl")
    include("plots.jl")
end
