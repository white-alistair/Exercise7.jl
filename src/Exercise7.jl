module Exercise7 
    export SIR, SIRV, SIRVSeasonal, SIRVDecayingImmunity, plot

    using DynamicalSystems, CairoMakie, Printf

    include("models.jl")
    include("plots.jl")
end
