using Exercise7
using Test

@testset verbose = true "Exercise7.jl" begin
    include("SIR.jl")
    include("SIRV.jl")
end
