using Test
using TropicalTensors

@testset "tropical tensor" begin
    function get_TroI()
        I2 = zeros(2,2)
        I2 .= -Inf
        I2[1,1] = 0.0
        I2[2,2] = 0.0

        I3 = zeros(2,2,2)
        I3 .= -Inf
        I3[1,1,1] = 0.0
        I3[2,2,2] = 0.0

        I4 = zeros(2,2,2,2)
        I4 .= -Inf
        I4[1,1,1,1] = 0.0
        I4[2,2,2,2] = 0.0

        return Tropical.(I2),Tropical.(I3),Tropical.(I4)
    end
    I2, I3, I4 = get_TroI()
    @test hypercubicI(Tropical{Float64}, 2, 2) ≈ I2
    @test hypercubicI(Tropical{Float64}, 3, 2) ≈ I3
    @test hypercubicI(Tropical{Float64}, 4, 2) ≈ I4
    @test spinglass_vertex_tensor(4) ≈ I4

    _get_J(::Val{:ferro}) = 1.0
    function get_TroB(v::Val)
        Jij = _get_J(v)
        tij = Tropical(Jij)
        _tij = Tropical(-Jij)
        return [tij _tij; _tij tij]
    end
    @test spinglass_bond_tensor(1.0) ≈ get_TroB(Val{:ferro}())
end
