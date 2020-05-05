using TropicalTensors
using LinearAlgebra, Test

@testset "i-bond tensor" begin
    mat = ones(Tropical{Float64}, 2, 2)
    J = 0.3
    @test spinglass_bond_tensor!(mat, J)[1] ≈ spinglass_bond_tensor(J)
    mat = Diagonal(ones(Tropical{Float64}, 4))
    @test spinglass_g4_tensor!(mat, J)[1] ≈ spinglass_g4_tensor(J)
    mat = ones(Tropical{Float64}, 16, 16)
    Js = randn(16)
    @test spinglass_g16_tensor!(mat, Js)[1] ≈ spinglass_g16_tensor(Js)
end
