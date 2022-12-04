function plot(model::SIR{T}; endtime::T = 365.0, Δt::T = 1.0, u0 = rand(T, 3)) where {T<:AbstractFloat}
    u0 /= sum(u0)  # Normalise so that S + I + R = 1
    ds = ContinuousDynamicalSystem(model, u0, nothing)
    tr = trajectory(ds, endtime; Δt)
    S, I, R = columns(tr)

    @printf "Initial conditions: S(0) = %.2f, I(0) = %.2f, R(0) = %.2f\n" u0[1] u0[2] u0[3]
    @printf "Final state:        S(T) = %.2f, I(T) = %.2f, R(T) = %.2f" S[end] I[end] R[end]

    # Plotting
    fig = Figure(resolution=(1200, 600))
    ax = Axis(fig[1, 1], xlabel="days", title= @sprintf "SIR Model with β = %.2f, γ = %.3f" model.β model.γ)
    times = Δt * collect(1:length(S))
    lines!(ax, times, S, linewidth=3, label="S")
    lines!(ax, times, I, linewidth=3, label="I")
    lines!(ax, times, R, linewidth=3, label="R")
    axislegend()
    fig
end

function plot(model::SIRV{T}; endtime::T = 100.0, Δt::T = 1.0, u0 = rand(T, 4)) where {T<:AbstractFloat}
    u0 /= sum(u0)  # Normalise so that S + I + R + V = 1
    ds = ContinuousDynamicalSystem(model, u0, nothing)
    tr = trajectory(ds, endtime; Δt)
    S, I, R, V = columns(tr)

    @printf "Initial conditions: S(0) = %4i, I(0) = %4i, R(0) = %4i, V(0) = %4i\n" u0[1] u0[2] u0[3] u0[4]
    @printf "Final state:        S(T) = %4i, I(T) = %4i, R(T) = %4i, V(T) = %4i" S[end] I[end] R[end] V[end]

    # Plotting
    fig = Figure(resolution=(1200, 600))
    ax = Axis(fig[1, 1], xlabel="days", title= @sprintf "SIRV Model with β = %.2f, γ = %.3f, ν = %.2f" model.β model.γ model.ν)
    times = Δt * collect(1:length(S))
    lines!(ax, times, S, linewidth=3, label="S")
    lines!(ax, times, I, linewidth=3, label="I")
    lines!(ax, times, R, linewidth=3, label="R")
    lines!(ax, times, V, linewidth=3, label="V")
    axislegend()
    fig
end
