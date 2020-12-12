using Test
using TropicalTensors
using CUDA, CuYao
using ForwardDiff
using Random

CUDA.allowscalar(false)

@testset "cuda" begin
    T = Tropical{Float64}
    @test ones(T, 3) == Array(CUDA.ones(T, 3))
    @test zeros(T, 3) == Array(CUDA.zeros(T, 3))
end

@testset "spinglass" begin
    lt = SquareLattice(10, 10)
    sg = Spinglass(lt, ones(Float32, 180), zeros(Float32, 100))
    res = solve(sg; usecuda=true)
    @test res.n == 180
end

@testset "test Chimera" begin
    lt = ChimeraLattice(3, 3)
    sg = Spinglass(lt, ones(Float32, 12*4 + 9*16), zeros(Float32, 9*8))
    res = solve(sg; usecuda=true)
    @test res.n == 12*4 + 9*16
end

@testset "forwarddiff-gpu" begin
    L = 5
    Js = ones(Float32, 2L*(L-1))
    gs = ForwardDiff.gradient(x->solve(SquareLattice(L, L), x, zeros(eltype(x), L^2); usecuda=true).n, Js)
    @test gs ≈ ones(Float32, 2L*(L-1))
end

@testset "counting tropical matmul" begin
    a = CountingTropical{Float64}.(randn(5, 7))
    b = CountingTropical{Float64}.(randn(7, 5))
    @test all(a * b .≈ Array(CuArray(a) * CuArray(b)))
end
