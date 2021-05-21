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
    for i=1:10
        lt = SquareLattice(10, 10)
        sg = Spinglass(lt, randn(Float32, 180), zeros(Float32, 100))
        res = solve(sg; usecuda=true)
        res0 = solve(sg; usecuda=false)
        @test res.n == res0.n
    end
end

@testset "masked square" begin
    for i=1:10
        lt = rand_maskedsquare(5, 5, 0.7)
        sg = Spinglass(lt, randn(Float32, length(sgbonds(lt))), zeros(Float32, length(lt)))
        res = solve(sg; usecuda=true)
        res0 = solve(sg; usecuda=false)
        @test res.n == res0.n
    end
end

@testset "test Chimera" begin
    for i=1:10
        lt = ChimeraLattice(3, 3)
        sg = Spinglass(lt, randn(Float32, 12*4 + 9*16), zeros(Float32, 9*8))
        res = solve(sg; usecuda=true)
        res0 = solve(sg; usecuda=false)
        @test res.n == res0.n
    end
end

@testset "forwarddiff-gpu" begin
    L = 5
    Js = ones(Float32, 2L*(L-1))
    gs = ForwardDiff.gradient(x->solve(SquareLattice(L, L), x, zeros(eltype(x), L^2); usecuda=true).n, Js)
    @test gs ≈ ones(Float32, 2L*(L-1))
end

@testset "counting tropical matmul" begin
    a = CountingTropicalF64.(randn(5, 7))
    b = CountingTropicalF64.(randn(7, 5))
    @test all(a * b .≈ Array(CuArray(a) * CuArray(b)))
end
