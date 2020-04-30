using Test
using TropicalTensors

@testset "tropical" begin
    @test Tropical(3) * Tropical(4) == Tropical(7)
    @test Tropical(3) + Tropical(4) == Tropical(4)
    @test Tropical(4) + Tropical(-1) == Tropical(4)
    @test zero(Tropical(2)) == Tropical(-999999)
    @test zero(Tropical(2.0)) == Tropical(-Inf)
    @test one(Tropical(2)) == Tropical(0)
    @test Tropical(2.0) ≈ Tropical(2.0 + 1e-10)
    @test Tropical(2) ≈ Tropical(2.0)
end

using CuArrays
@testset "cuda" begin
    T = Tropical{Float64} 
    @test ones(T, 3) == Array(CuArrays.ones(T, 3))
    @test zeros(T, 3) == Array(CuArrays.zeros(T, 3))
end 
