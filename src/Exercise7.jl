module Exercise7 
    export SIR, SIRV, SIRVSeasonal, SIRVDecayingImmunity, plot, plot_vaccination_rate

    using DynamicalSystems, CairoMakie, Printf

    include("models.jl")
    include("plots.jl")
end
