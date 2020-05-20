using TropicalTensors
using Yao, TropicalYao
using Test
using Viznet
using NiLang, NiLang.AD
using Random

@testset "spinglass" begin
    Random.seed!(2)
    lt = rand_maskedsquare(8, 10, 0.8)
    sg = Spinglass(lt, rand(Int32[-1,1], length(sgbonds(lt))), zeros(Int32, length(lt)))
    Lx, Ly = size(sg.lattice)
    reg = ArrayReg(ones(Tropical{Int32}, 1<<(Ly+2)))
    STK = stack4reg(reg, 100)
    @test isolve_largemem(Int32(0), sg, reg, STK)[1] == solve(sg).n
    A = stack4reg(reg, 2Ly)
    B = stack4reg(reg, Lx-1)
    @test check_inv(isolve, (Int32(0), sg, reg, A, B))
    reg = ArrayReg(ones(Tropical{Int32}, 1<<(Ly+2)))
    A = stack4reg(reg, 2Ly)
    B = stack4reg(reg, Lx-1)
    @test isolve(Int32(0), sg, reg, A, B)[1] == solve(sg).n
end

@testset "optconfig" begin
    Random.seed!(5)
    function opt_config_largemem(sg::Spinglass{LT,T}) where {LT,T}
        reg = TropicalTensors._init_reg(T, sg.lattice, Val(:false))
        A = stack4reg(reg, length(sg.lattice)*2)
        eng, sg, reg, A = isolve_largemem(T(0.0), sg, reg, A)
        sgg = Spinglass(sg.lattice, GVar.(sg.Js, zero(sg.Js)), GVar.(sg.hs, zero(sg.hs)))
        gres = (~isolve_largemem)(GVar(eng, T(1)), sgg, GVar(reg), GVar(A))
        return TropicalTensors.SpinglassOptConfig(sg, eng, grad.(gres[2].Js), grad.(gres[2].hs))
    end
    sg = rand_spinglass(Int32, rand_maskedsquare(9, 7, 0.8); jt=Randpm(), ht=Zero(), seed=2)
    optc = opt_config(sg)
    optc2 = opt_config_largemem(sg)
    @test optc.grad_Js == optc2.grad_Js
    @test optc.grad_hs == optc2.grad_hs
    @test sum(sg.Js .* optc.grad_Js) == optc.eng
    @test sum(sg.Js .* optc2.grad_Js) == optc2.eng
    @test optc.eng == optc2.eng
end

@testset "fallback" begin
    function opt_config_largemem(sg::Spinglass{LT,T}) where {LT,T}
        reg = TropicalTensors._init_reg(T, sg.lattice, Val(:false))
        A = stack4reg(reg, length(sg.lattice)*2)
        eng, sg, reg, A = isolve_largemem(T(0.0), sg, reg, A)
        sgg = Spinglass(sg.lattice, GVar.(sg.Js, zero(sg.Js)), GVar.(sg.hs, zero(sg.hs)))
        gres = (~isolve_largemem)(GVar(eng, T(1)), sgg, GVar(reg), GVar(A))
        return TropicalTensors.SpinglassOptConfig(sg, eng, grad.(gres[2].Js), grad.(gres[2].hs))
    end
    sg = Spinglass(MaskedSquareLattice(ones(Bool,2,2)), [-1,1,1,1,1,-1], [0,0,0,0])
    Spinglass{MaskedSquareLattice,Int64}(MaskedSquareLattice(Bool[1 1; 1 1]), [-1, 1, 1, 1, 1, -1], [0, 0, 0, 0])
    res = opt_config_largemem(sg)
    @test res.eng == 2
end
