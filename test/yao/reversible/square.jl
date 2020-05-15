using TropicalTensors
using Yao, TropicalYao
using Test
using Viznet
using NiLang, NiLang.AD

@testset "spinglass" begin
    sg = Spinglass(SquareLattice(8, 10), ones(Int32, 142), zeros(Int32, 80))
    Lx, Ly = size(sg.lattice)
    reg = ArrayReg(ones(Tropical{Int32}, 1<<Ly))
    A = stack4reg(reg, Ly)
    B = stack4reg(reg, Lx-1)
    @test check_inv(isolve, (Int32(0), sg, reg, A, B))
    reg = ArrayReg(ones(Tropical{Int32}, 1<<Ly))
    A = stack4reg(reg, Ly)
    B = stack4reg(reg, Lx-1)
    @test isolve(Int32(0), sg, reg, A, B)[1] === Int32(142)
end

@testset "optconfig" begin
    function opt_config_largemem(sg::Spinglass{LT,T}) where {LT,T}
        reg = ArrayReg(ones(Tropical{T}, 1<<sg.lattice.Ny))
        A = stack4reg(reg, length(sg.lattice))
        eng, sg, reg, A = isolve_largemem(T(0.0), sg, reg, A)
        sgg = Spinglass(sg.lattice, GVar.(sg.Js, zero(sg.Js)), GVar.(sg.hs, zero(sg.hs)))
        gres = (~isolve_largemem)(GVar(eng, T(1)), sgg, GVar(reg), GVar(A))
        empty!(NiLang.GLOBAL_STACK)
        return TropicalTensors.SpinglassOptConfig(sg, eng, grad.(gres[2].Js), grad.(gres[2].hs))
    end
    sg = rand_spinglass(Int32, SquareLattice(9, 7); jt=Randpm(), ht=Zero(), seed=2)
    optc = opt_config(sg)
    optc2 = opt_config_largemem(sg)
    @test optc.grad_J == optc2.grad_J
    @test optc.eng == optc2.eng
end
