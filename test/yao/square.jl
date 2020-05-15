using TropicalTensors
using Test

@testset "bonds" begin
    lt = SquareLattice(3,2)
    @test sgbonds(lt) == [(1, 4), (1, 2), (4, 5), (2, 5), (2, 3), (5, 6), (3, 6)]
end

@testset "spinglass" begin
    lt = SquareLattice(10, 8)
    sg = Spinglass(lt, ones(Float32, 142), zeros(Float32, 80))
    res = solve(sg; usecuda=false)
    @test res.n == 142
end
