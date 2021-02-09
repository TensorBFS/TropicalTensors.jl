using TropicalTensors
using Test, Random
using Viznet
using Yao
using TropicalYao: hypercubicI

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
    @test_broken p1 ≈ p2
    # sum-reset gate
    function greset(::Type{T}) where T
        matblock([one(T) one(T); zero(T) zero(T)])
    end
    reg |> put(8, 7=>greset(Float64))
    @test measure(reg, 7; nshots=10) == zeros(10)
    focus!(reg, (2,))
    p3 = probs(reg)
    relax!(reg, (2,))
    @test_broken p3 ≈ p1
end

@testset "MaskedSquareLattice" begin
    lt = rand_maskedsquare(4, 4, 1.0)
    sq = SquareLattice(4, 4)
    @test length(TropicalTensors.sgvertices(lt)) == 16
    @test length(TropicalTensors.sgbonds(lt)) == length(bonds(sq, nth=1)) + length(bonds(sq, nth=2))
end

@testset "solve" begin
    lt = rand_maskedsquare(4, 4, 1.0)
    sg = rand_spinglass(Float32, lt; jt=Ferro(), ht=Zero())
    res = solve(sg)
    @test res.n == length(sgbonds(lt))
end

@testset "solve masked" begin
    lt = MaskedSquareLattice([1 1 1 1;
        1 1 0 1;
        1 1 1 1;
        1 1 1 0])
    sg = rand_spinglass(Float64, lt; jt=Ferro(), ht=Zero())
    res = solve(sg)
    @test res.n == 31
end

@testset "fallback cut" begin
    lt = MaskedSquareLattice([0 1; 1 0; 0 1])
    sg = Spinglass(lt, [-1,1], [0,0,0])
    @test solve(sg).n == 2
end

@testset "forwarddiff" begin
    L = 6
    lt = rand_maskedsquare(L, L, 0.8)
    hs = zeros(length(lt))
    Js = rand([-1,1.0], length(sgbonds(lt)))
    res = solve(lt, Js, hs; usecuda=false).n
    gs = ForwardDiff.gradient(x->(dn = solve(lt, eltype(x).(Js), x; usecuda=false).n; (@test ForwardDiff.value(dn)==res); dn), hs)
    function compute_loss(lt, Js, config)
        res = 0.0
        vs = vertices(lt)
        for (i, (src, dst)) in enumerate(sgbonds(lt))
            res += Js[i] * config[findfirst(==(src), vs)] * config[findfirst(==(dst), vs)]
        end
        return res
    end
    @test compute_loss(lt, Js, gs) == res
end
