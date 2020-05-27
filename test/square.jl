using TropicalTensors
using ForwardDiff
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

@testset "forwarddiff" begin
    L = 3
    Js = ones(Float32, 2L*(L-1))
    gs = ForwardDiff.gradient(x->solve(SquareLattice(L, L), x, zeros(eltype(x), L^2); usecuda=false).n, Js)
    @test gs ≈ ones(Float32, 2L*(L-1))
end
