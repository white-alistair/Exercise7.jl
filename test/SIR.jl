@testset "SIR.jl" begin
    @testset "Standard SIR Model" begin
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
end
