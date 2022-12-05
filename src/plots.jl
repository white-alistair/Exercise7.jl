function plot(
    model::AbstractSIRModel{T};
    endtime::T = 365.0,
    Δt::T = 1.0,
    u0::Vector{T} = rand(T, 3),
) where {T<:AbstractFloat}
    u0 /= sum(u0)  # Normalise so that S + I + R = 1
    ds = ContinuousDynamicalSystem(model, u0, nothing)
    tr = trajectory(ds, endtime; Δt)
    S, I, R = columns(tr)

    @printf "Initial conditions: S(0) = %.2f, I(0) = %.2f, R(0) = %.2f\n" u0[1] u0[2] u0[3]
    @printf "Final state:        S(T) = %.2f, I(T) = %.2f, R(T) = %.2f" S[end] I[end] R[end]

    # Plotting
    fig = Figure(; resolution = (1200, 600))
    ax = Axis(fig[1, 1]; xlabel = "days", title = repr(model))
    times = Δt*1:length(S)
    lines!(ax, times, S; linewidth = 3, label = "S")
    lines!(ax, times, I; linewidth = 3, label = "I")
    lines!(ax, times, R; linewidth = 3, label = "R")
    axislegend(ax)
    return fig
end

function plot(
    model::AbstractSIRVModel{T};
    endtime::T = 365.0,
    Δt::T = 1.0,
    u0::Vector{T} = rand(T, 4),
) where {T<:AbstractFloat}
    u0 /= sum(u0)  # Normalise so that S + I + R + V = 1
    ds = ContinuousDynamicalSystem(model, u0, nothing)
    tr = trajectory(ds, endtime; Δt)
    S, I, R, V = columns(tr)

    @printf "Initial conditions: S(0) = %.2f, I(0) = %.2f, R(0) = %.2f, V(0) = %.2f\n" u0[1] u0[2] u0[3] u0[4]
    @printf "Final state:        S(T) = %.2f, I(T) = %.2f, R(T) = %.2f, V(T) = %.2f" S[end] I[end] R[end] V[end]

    # Plotting
    fig = Figure(; resolution = (1200, 600))
    ax = Axis(fig[1, 1]; xlabel = "days", title = repr(model))
    times = Δt*1:length(S)
    lines!(ax, times, S; linewidth = 3, label = "S")
    lines!(ax, times, I; linewidth = 3, label = "I")
    lines!(ax, times, R; linewidth = 3, label = "R")
    lines!(ax, times, V; linewidth = 3, label = "V")
    axislegend(ax)
    return fig
end

function plot_phase_diagram(
    model::AbstractSIRVModel{T};
    endtime::T = 1000.0,
    Ttr::T = 1000.0,
    u0::Vector{T} = rand(4),
    tolerance::T = 1e-15,
) where {T}
    # Initial conditions
    u0 /= sum(u0)  # Normalise so that S + I + R + V = 1
    @printf "Initial conditions: S(0) = %.2f, I(0) = %.2f, R(0) = %.2f, V(0) = %.2f" u0[1] u0[2] u0[3] u0[4]

    # Integration parameters
    cb = PositiveDomain(zeros(T, 4); abstol = eps(T))  # Callback to ensure the solution remains positive

    # Transient integration
    prob_transient = ODEProblem(model, u0, (zero(T), Ttr))
    sol_transient = solve(
        prob_transient,
        Vern9();
        abstol = tolerance,
        reltol = tolerance,
        callback = cb,
    )

    # Actual integration
    u1 = sol_transient[end]
    prob = ODEProblem(model, u1, (zero(T), endtime))
    sol = solve(
        prob,
        Vern9();
        saveat = 1.0,
        abstol = tolerance,
        reltol = tolerance,
        callback = cb,
    )

    # Plotting
    S = sol[1, :]
    I = sol[2, :]
    fig = Figure(; resolution = (1200, 600))
    ax1 = Axis(fig[1, 1]; xlabel = "day", ylabel = "-log(I)")
    ax2 = Axis(fig[1, 2]; xlabel = "-log(S)", ylabel = "-log(I)")
    lines!(ax1, sol.t, -1 * log.(I); linewidth = 3)
    lines!(ax2, -1 * log.(S), -1 * log.(I); linewidth = 3)
    return fig
end

function plot_total_infections_by_vax_rate(
    ν_range::AbstractVector{T};
    endtime::T = 365.0,
    Δt::T = 1.0,
    u0::Vector{T} = rand(T, 4),
) where {T<:AbstractFloat}
    total_infections = T[]

    for ν in ν_range
        model = SIRV(ν)
        ds = ContinuousDynamicalSystem(model, u0, nothing)
        tr = trajectory(ds, endtime; Δt)
        R_end = tr[end][3]
        push!(total_infections, R_end)
    end

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xlabel = "Vaccination Rate",
        xticks = ν_range[1]:0.1:ν_range[end],
        xminorticks = IntervalsBetween(10),
        xminorgridvisible = true,
        ylabel = "Total Infections Before Extinction",
        title = "SIRV Model with u0 = $u0",
    )
    scatter!(ax, ν_range, total_infections)
    return fig
end
