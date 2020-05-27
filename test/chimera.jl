using TropicalTensors
using Test
using Yao

@testset "test red reg" begin
    reg = TropicalTensors.red_reg(Float32, 3, randn(8+16*3), randn(4*3); usecuda=false)
    @test reg isa ArrayReg{1,Tropical{Float32},Matrix{Tropical{Float32}}}
end

@testset "sg bonds" begin
    lt = ChimeraLattice(3, 2)
    bonds = sgbonds(lt)
    bonds[1] == (1, 25)
    bonds[2] == (1, 26)
    bonds[5] == (5, 1)
    bonds[37] == (5, 13)
    @test length(bonds) == 16*6 + 3*4 + 16
end

@testset "test Chimera" begin
    lt = ChimeraLattice(3, 3)
    sg = Spinglass(lt, ones(Float32, 12*4 + 9*16), zeros(Float32, 9*8))
    res = solve(sg; usecuda=false)
    @test res.n == 12*4 + 9*16
end
