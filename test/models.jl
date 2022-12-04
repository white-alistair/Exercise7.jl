@testset "models.jl" begin
    @testset "SIR Model" begin
        # Constructors
        β, γ = 1.0, 2.0
        model = SIR(β, γ)
        @test model.β == β
        @test model.γ == γ

        # Defaults
        model = SIR()
        @test model.β == 0.35
        @test model.γ == 0.035

        # Equations
        du = zeros(3)
        u = [1.0, 2.0, 3.0]
        model(du, u, nothing, nothing)
        @test du ≈ [-0.7, 0.63, 0.07]
    end

    @testset "SIR Model with Decaying Immunity" begin
        # Constructors
        β, γ, α = 1.0, 2.0, 3.0
        model = SIRDecayingImmunity(β, γ, α)
        @test model.β == β
        @test model.γ == γ
        @test model.α == α

        # Defaults
        model = SIRDecayingImmunity(α)
        @test model.β == 0.35
        @test model.γ == 0.035
        @test model.α == α

        # Equations
        du = zeros(3)
        u = [1.0, 2.0, 3.0]
        model(du, u, nothing, nothing)
        @test du ≈ [8.3, 0.63, -8.93]
    end

    @testset "SIRV Model" begin
        # Constructors
        β, γ, ν = 1.0, 2.0, 3.0
        model = SIRV(β, γ, ν)
        @test model.β == β
        @test model.γ == γ
        @test model.ν == ν

        # Defaults
        model = SIRV(ν)
        @test model.β == 0.35
        @test model.γ == 0.035
        @test model.ν == ν

        # Equations
        du = zeros(4)
        u = [1.0, 2.0, 3.0, 4.0]
        model(du, u, nothing, nothing)
        @test du ≈ [-3.7, 0.63, 0.07, 3.0]
    end

    @testset "SIRV Model with Seasonal Contact" begin
        # Constructors
        β₀, β₁, γ, ν = 1.0, 2.0, 3.0, 4.0
        model = SIRVSeasonalContact(β₀, β₁, γ, ν)
        @test model.β₀ == β₀
        @test model.β₁ == β₁
        @test model.γ == γ
        @test model.ν == ν

        # Defaults
        model = SIRVSeasonalContact(1.0, 2.0)
        @test model.β₀ == 0.35
        @test model.β₁ == 1.0
        @test model.γ == 0.035
        @test model.ν == 2.0

        # Equations
        du = zeros(4)
        u = [1.0, 2.0, 3.0, 4.0]
        model(du, u, nothing, 0.0)
        @test du ≈ [-3.4, 1.33, 0.07, 2.0]
    end

    @testset "SIRV Model with Decaying Immunity" begin
        # Constructors
        β, γ, ν, α, μ = 1.0, 2.0, 3.0, 4.0, 5.0
        model = SIRVDecayingImmunity(β, γ, ν, α, μ)
        @test model.β == β
        @test model.γ == γ
        @test model.ν == ν
        @test model.α == α
        @test model.μ == μ

        # Defaults
        model = SIRVDecayingImmunity(ν, α, μ)
        @test model.β == 0.35
        @test model.γ == 0.035
        @test model.ν == ν
        @test model.α == α
        @test model.μ == μ

        # Equations
        du = zeros(4)
        u = [1.0, 2.0, 3.0, 4.0]
        model(du, u, nothing, nothing)
        @test du ≈ [28.3, 0.63, -11.93, -17.0]
    end

    @testset "SIRV Model with Seasonal Contact and Decaying Immunity" begin
        # Constructors
        β₀, β₁, γ, ν, α, μ = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0
        model = SIRVSeasonalContactDecayingImmunity(β₀, β₁, γ, ν, α, μ)
        @test model.β₀ == β₀
        @test model.β₁ == β₁
        @test model.γ == γ
        @test model.ν == ν
        @test model.α == α
        @test model.μ == μ

        # Defaults
        model = SIRVSeasonalContactDecayingImmunity(β₁, ν, α, μ)
        @test model.β₀ == 0.35
        @test model.β₁ == β₁
        @test model.γ == 0.035
        @test model.ν == ν
        @test model.α == α
        @test model.μ == μ

        # Equations
        du = zeros(4)
        u = [1.0, 2.0, 3.0, 4.0]
        model(du, u, nothing, 0.0)
        @test du ≈ [32.9, 2.03, -14.93, -20.0]
    end
end
