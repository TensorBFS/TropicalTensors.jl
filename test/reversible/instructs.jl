using TropicalTensors.Reversible, TropicalTensors
using Test
using LinearAlgebra
using Yao

@testset "new instr" begin
    g4 = Diagonal(spinglass_g4_tensor(1.5))
    reg = ArrayReg(randn(1<<12) .|> Tropical)
    s1 = i_instruct!(copy(vec(reg.state)), g4, (7, 3), (), ())[1]
    nreg = copy(reg) |> put(12, (7, 3)=>matblock(g4))
    @test statevec(nreg) ≈ s1

    g4 = randn(4, 4) .|> Tropical
    reg = ArrayReg(randn(1<<12) .|> Tropical)
    s1 = i_instruct!(copy(vec(reg.state)), g4, (7, 3), (), ())[1]
    nreg = copy(reg) |> put(12, (7, 3)=>matblock(g4))
    @test statevec(nreg) ≈ s1
end
