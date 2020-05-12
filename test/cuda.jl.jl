using Test
using TropicalTensors

using CuArrays
@testset "cuda" begin
    T = Tropical{Float64}
    @test ones(T, 3) == Array(CuArrays.ones(T, 3))
    @test zeros(T, 3) == Array(CuArrays.zeros(T, 3))
end
