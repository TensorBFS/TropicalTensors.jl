using TropicalTensors, OMEinsum
using Test, Random
using Viznet
using Yao

@testset "copy vertex" begin
    m = hypercubicI(3, 2)
    ct = copyvertex(Float64)
    cm = reshape(ein"kal,iaj->ikjl"(ct, m), 4, 4)
    cm_adjoint = reshape(ein"kal,iaj->ikjl"(permutedims(copyvertex(Float64), (3,2,1)), m), 4, 4)
    @test cm == copytensor(Float64)
    @test cm_adjoint == cm'
end

@testset "copy and reset gates" begin
    reg = join(zero_state(Float64, 2), rand_state(Float64, 6))
    focus!(reg, (2,))
    p1 = probs(reg)
    relax!(reg, (2,))
    @test measure(reg, 7; nshots=10) == zeros(10)

    # copy information from 2 to 7
    function gcp(::Type{T}) where T
        matblock(copytensor(T))
    end
    reg |> put(8, (2,7)=>gcp(Float64))
    focus!(reg, (7,))
    p2 = probs(reg)
    relax!(reg, (7,))
    @test p1 ≈ p2
    # sum-reset gate
    function greset(::Type{T}) where T
        matblock([one(T) one(T); zero(T) zero(T)])
    end
    reg |> put(8, 7=>greset(Float64))
    @test measure(reg, 7; nshots=10) == zeros(10)
    focus!(reg, (2,))
    p3 = probs(reg)
    relax!(reg, (2,))
    @test p3 ≈ p1
end

@testset "SquareLattice2nd" begin
    lt = SquareLattice2nd(4, 4)
    sq = SquareLattice(4, 4)
    @test length(TropicalTensors.sgvertexorder(lt)) == 16
    @test length(TropicalTensors.sgbonds(lt)) == length(bonds(sq, nth=1)) + length(bonds(sq, nth=2))
end

@testset "solve" begin
    lt = SquareLattice2nd(4, 4)
    sg = rand_spinglass(Float64, lt; jt=Ferro(), ht=Zero())
    res = solve(sg)
    @test res.n == length(sgbonds(lt))
end
