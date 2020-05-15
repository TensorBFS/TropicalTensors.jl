using TropicalTensors, TropicalTensors.Reversible
using Test

@testset "spinglass" begin
    L = 10
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    A = stack4reg(reg, L)
    B = stack4reg(reg, L-1)
    @test check_inv(spinglass_yao, (Float32(0.0), reg, L, ones(Float32, 180), A, B))
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    A = stack4reg(reg, L)
    B = stack4reg(reg, L-1)
    @test spinglass_yao(Float32(0.0), reg, L, ones(Float32,180), A, B)[1] ≈ 180.0
    empty!(NiLang.GLOBAL_STACK)
end
