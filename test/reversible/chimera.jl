using TropicalTensors, TropicalTensors.Reversible
using Yao, TropicalYao
using Test
using Viznet
using NiLang, NiLang.AD

@testset "isolve" begin
    Lx, Ly = 3, 4
    lt = ChimeraLattice(Lx, Ly)
    nJ = length(bonds(lt))
    sg = Spinglass(lt, ones(Int32, nJ), zeros(Int32, Lx*Ly*8))
    reg = ArrayReg(ones(Tropical{Int32}, 1<<(4*Ly)))
    A = stack4reg(reg, 4*Ly)
    B = stack4reg(reg, Lx-1)
    @test check_inv(isolve, (Int32(0), sg, reg, A, B))
    reg = ArrayReg(ones(Tropical{Int32}, 1<<(4*Ly)))
    A = stack4reg(reg, 4*Ly)
    B = stack4reg(reg, Lx-1)
    @test isolve(Int32(0), sg, reg, A, B)[1] === Int32(nJ)
end

function _solve(Lx, Ly; largemem)
    lt = ChimeraLattice(Lx, Ly)
    nJ = length(bonds(lt))
    sg = Spinglass(lt, ones(Int32, nJ), zeros(Int32, Lx*Ly*8))
    reg = ArrayReg(ones(Tropical{Int32}, 1<<(4*Ly)))
    if largemem
        A = stack4reg(reg, Ly*Lx + Ly*(Lx-1)*4 + Lx-1)
        isolve_largemem(Int32(0), sg, reg, A)[1]
    else
        A = stack4reg(reg, Ly*4)
        B = stack4reg(reg, Lx-1)
        isolve(Int32(0), sg, reg, A, B)[1]
    end
end

@testset "isolve largemem" begin
    lt = ChimeraLattice(3, 2)
    _, Lx, Ly = size(lt)
    nJ = length(bonds(lt))
    sg = Spinglass(lt, ones(Int32, nJ), zeros(Int32, Lx*Ly*8))
    reg = ArrayReg(ones(Tropical{Int32}, 1<<(4*Ly)))
    A = stack4reg(reg, Ly*Lx + Ly*(Lx-1)*4 + Lx-1)
    @test check_inv(isolve_largemem, (Int32(0), sg, reg, A))
    reg = ArrayReg(ones(Tropical{Int32}, 1<<(4*Ly)))
    A = stack4reg(reg, Ly*Lx + Ly*(Lx-1)*4 + Lx-1)
    @test isolve_largemem(Int32(0), sg, reg, A)[1] === Int32(nJ)
end

@testset "optconfig" begin
    function opt_config_largemem(sg::Spinglass{LT,T}) where {LT,T}
        Ly, Lx = sg.lattice.Ny, sg.lattice.Nx
        reg = ArrayReg(ones(Tropical{T}, 1<<(Ly*4)))
        A = stack4reg(reg, Ly*Lx + Ly*(Lx-1)*4 + Lx-1)
        eng, sg, reg, A = isolve_largemem(T(0.0), sg, reg, A)
        sgg = Spinglass(sg.lattice, GVar.(sg.Js, zero(sg.Js)), GVar.(sg.hs, zero(sg.hs)))
        gres = (~isolve_largemem)(GVar(eng, T(1)), sgg, GVar(reg), GVar(A))
        return TropicalTensors.SpinglassOptConfig(sg, eng, grad.(gres[2].Js), grad.(gres[2].hs))
    end
    sg = rand_spinglass(Int32, ChimeraLattice(3, 4); jt=Randpm(), ht=Zero(), seed=2)
    optc = opt_config(sg)
    optc2 = opt_config_largemem(sg)
    @test optc.grad_Js == optc2.grad_Js
    @test optc.grad_hs == optc2.grad_hs
    @test optc.eng == optc2.eng
end
