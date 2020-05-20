using TropicalTensors
using Yao, TropicalYao
using Test
using Viznet
using NiLang, NiLang.AD
using Random
using TropicalTensors: cachesize_A, cachesize_B, cachesize_largemem

@testset "spinglass" begin
    Random.seed!(2)
    lt = rand_maskedsquare(8, 10, 0.8)
    sg = Spinglass(lt, rand(Int32[-1,1], length(sgbonds(lt))), zeros(Int32, length(lt)))
    Lx, Ly = size(sg.lattice)
    reg = ArrayReg(ones(Tropical{Int32}, 1<<(Ly+2)))
    STK = stack4reg(reg, cachesize_largemem(lt))
    @test isolve_largemem(Int32(0), sg, reg, STK)[1] == solve(sg).n
    A = stack4reg(reg, cachesize_A(lt))
    B = stack4reg(reg, cachesize_B(lt))
    @test check_inv(isolve, (Int32(0), sg, reg, A, B))
    reg = ArrayReg(ones(Tropical{Int32}, 1<<(Ly+2)))
    A = stack4reg(reg, cachesize_A(lt))
    B = stack4reg(reg, cachesize_B(lt))
    @test isolve(Int32(0), sg, reg, A, B)[1] == solve(sg).n
end

@testset "optconfig" begin
    Random.seed!(5)
    sg = rand_spinglass(Int64, rand_maskedsquare(9, 7, 0.8); jt=Randpm(), ht=Zero(), seed=2)
    optc = opt_config(sg)
    optc2 = opt_config_largemem(sg)
    @test optc.grad_Js == optc2.grad_Js
    @test optc.grad_hs == optc2.grad_hs
    @test sum(sg.Js .* optc.grad_Js) == optc.eng
    @test sum(sg.Js .* optc2.grad_Js) == optc2.eng
    @test optc.eng == optc2.eng
end

@testset "fallback" begin
    sg = Spinglass(MaskedSquareLattice(ones(Bool,2,2)), [-1,1,1,1,1,-1], [0,0,0,0])
    res = opt_config_largemem(sg)
    @test res.eng == 2
end
