using Test
using TropicalTensors

@testset "tropical" begin
    @test Tropical(3) * Tropical(4) == Tropical(7)
    @test Tropical(3) + Tropical(4) == Tropical(4)
    @test Tropical(4) + Tropical(-1) == Tropical(4)
    @test zero(Tropical(2)) == Tropical(-9223372036854775808)
    @test zero(Tropical(2.0)) == Tropical(-Inf)
    @test zero(Tropical{Float32}) == Tropical(-Inf32)
    @test zero(Tropical{Float64}) == Tropical(-Inf)
    @test zero(Tropical{Int16}) == Tropical(Int16(-32768))
    @test zero(Tropical{Int32}) == Tropical(Int32(-2147483648))
    @test zero(Tropical{Int64}) == Tropical(Int64(-9223372036854775808))
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
