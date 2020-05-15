using Test
using TropicalTensors
using CuArrays, CuYao
CuArrays.allowscalar(false)

@testset "cuda" begin
    T = Tropical{Float64}
    @test ones(T, 3) == Array(CuArrays.ones(T, 3))
    @test zeros(T, 3) == Array(CuArrays.zeros(T, 3))
end

@testset "spinglass" begin
    lt = SquareLattice(10, 10)
    sg = Spinglass(lt, ones(Float32, 180), zeros(Float32, 100))
    res = solve(sg; usecuda=true)
    @test res.n == 180
end

@testset "test Chimera" begin
    lt = ChimeraLattice(3, 3)
    sg = Spinglass(lt, ones(Float32, 12*4 + 9*16), zeros(Float32, 3*8))
    res = solve(sg; usecuda=true)
    @test res.n == 12*4 + 9*16
end
