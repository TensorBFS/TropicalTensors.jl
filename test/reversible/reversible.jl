using Test
using TropicalTensors, TropicalTensors.Reversible

@testset "reversible square" begin
    include("square.jl")
end

@testset "reversible chimera" begin
    include("chimera.jl")
end

@testset "reversible second_neighbor" begin
    include("second_neighbor.jl")
end
