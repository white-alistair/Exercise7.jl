# Exercise7.jl Documentation

## Compartmental Models
```@docs
Exercise7.AbstractCompartmentalModel
```

## SIR Models
```@docs
Exercise7.AbstractSIRModel
Exercise7.SIR
Exercise7.SIR(du, u, p, t)
Exercise7.SIR()

Exercise7.SIRDecayingImmunity
Exercise7.SIRDecayingImmunity(α)
Exercise7.SIRDecayingImmunity(du, u, p, t)
```

## SIRV Models
```@docs
Exercise7.AbstractSIRVModel

Exercise7.SIRV
Exercise7.SIRV(du, u, p, t)
Exercise7.SIRV(ν)

Exercise7.SIRVSeasonalContact
Exercise7.SIRVSeasonalContact(β₁, ν)
Exercise7.SIRVSeasonalContact(du, u, p, t)

Exercise7.SIRVDecayingImmunity
Exercise7.SIRVDecayingImmunity(ν, α, μ)
Exercise7.SIRVDecayingImmunity(du, u, p, t)

Exercise7.SIRVSeasonalContactDecayingImmunity
Exercise7.SIRVSeasonalContactDecayingImmunity(β₁, ν, α, μ)
```