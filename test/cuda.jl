using Test
using TropicalTensors
using CUDA, CuYao
using ForwardDiff
using TropicalTensors: c2l, l2c
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

@testset "c2l" begin
    for i=1:100
        shape = (4,rand(1:5),rand(1:7),5,19)
        target = ([rand(1:s) for s in shape]...,)
        @test c2l(shape, target) == LinearIndices(shape)[target...]
    end
    for i=1:100
        shape = (4,rand(1:5),rand(1:12),15,19)
        ci = CartesianIndices(shape)
        i = rand(1:prod(shape))
        @test l2c(shape, i) == ci[i].I
    end
end

@testset "permutedims" begin
    a = randn(rand(1:3, 20)...)
    A = CuArray(a)
    p = randperm(20)
    @test Array(permutedims(A, p)) ≈ permutedims(a, p)
end

@testset "counting tropical matmul" begin
    a = CountingTropical{Float64}.(randn(5, 7))
    b = CountingTropical{Float64}.(randn(7, 5))
    @test a * b ≈ Array(CuArray(a) * CuArray(b))
end
