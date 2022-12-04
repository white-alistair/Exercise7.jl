function plot(model::SIRV{T}; endtime::T = 1000.0, Δt::T = 0.1, u0 = 1000 * rand(T, 4)) where {T<:AbstractFloat}
    u0 /= sum(u0)  # Normalise so that S + I + R + V = 1
    ds = ContinuousDynamicalSystem(model, u0)
    tr = trajectory(ds, endtime; Δt)
    S, I, R, V = columns(tr)

    @printf "Initial conditions: S(0) = %4i, I(0) = %4i, R(0) = %4i, V(0) = %4i\n" u0[1] u0[2] u0[3] u0[4]
    @printf "Final state:        S(T) = %4i, I(T) = %4i, R(T) = %4i, V(T) = %4i" S[end] I[end] R[end] V[end]

    # Plotting
    fig = Figure(resolution=(1200, 600))
    ax = Axis(fig[1, 1], xlabel="days", title="SIRV Model")
    times = Δt * collect(1:length(S))
    lines!(ax, times, S, linewidth=3, label="S")
    lines!(ax, times, I, linewidth=3, label="I")
    lines!(ax, times, R, linewidth=3, label="R")
    lines!(ax, times, V, linewidth=3, label="V")
    axislegend()
    fig
end
