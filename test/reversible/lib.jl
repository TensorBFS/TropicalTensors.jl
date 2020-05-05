using TropicalTensors, TropicalTensors.Reversible
using Test, NiLang.AD, NiLang, Random

@testset "GVar" begin
    x = Tropical(0.4)
	@test (~Tropical)(x) == 0.4
    @test GVar(x) isa Tropical{<:GVar}
    @test grad(x) isa Tropical{Float64}
    @test (~GVar)(GVar(Tropical(0.4))) == x
end

@testset "Tropical algebra" begin
	Random.seed!(4)
	@test check_inv(tropical_muladd, (Tropical(2.0), Tropical(3.0), Tropical(4.0), false))
	@test check_inv(tropical_muladd, (Tropical(2.0), Tropical(3.0), Tropical(-5.0), false))
	x = Tropical.(randn(100, 100))
	y = Tropical.(randn(100, 100))
	v = Tropical.(randn(100))
	out1 = x * y
	out2 = one.(out1)
	igemm!(out2, x, y)
	out4 = Tropical.(zeros(100))
	out3 = x * v
	bk = zeros(Bool, 100)
	igemv!(out4, x, v, bk)
	@test out1 ≈ out2
	@test out3 ≈ out4
end

@testset "isum" begin
	isum(TropicalF32(0), TropicalF32.([3.0,10.0,2.0]))[1].n == 10
end
