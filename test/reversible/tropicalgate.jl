using TropicalTensors, TropicalTensors.Reversible
using Test, Random
using Yao
using NiLang

@testset "tropical block" begin
    reg = ArrayReg(Tropical.(randn(1<<10)))
    reg1 = apply_G2!(copy(reg), 2, 0.5)[1]
    reg2 = copy(reg) |> put(10, 2=>matblock(spinglass_bond_tensor(0.5)))
    @test reg1 ≈ reg2
    @test check_inv(apply_G2!, (copy(reg), 2, 0.5))
    reg1 = apply_G4!(copy(reg), (4,2), 0.5)[1]
    reg2 = copy(reg) |> put(10, (4,2)=>matblock(spinglass_g4_tensor(0.5)))
    @test reg1 ≈ reg2
    @test check_inv(apply_G4!, (copy(reg), (4,2), 0.5))
    Js = randn(16)
    reg1 = apply_G16!(copy(reg), (4,1,2,9), Js)[1]
    reg2 = copy(reg) |> put(10, (4,1,2,9)=>matblock(spinglass_g16_tensor(Js)))
    @test reg1 ≈ reg2
    @test check_inv(apply_G16!, (copy(reg), (4,1,2,9), Js))
end
