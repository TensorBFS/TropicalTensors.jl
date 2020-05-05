using TropicalTensors, TropicalTensors.Reversible
using Test

@testset "einsum" begin
	N = 3
	a = Tropical.(randn(N, N))
	b = Tropical.(randn(N, N))
	c = Tropical.(zeros(N, N))
	ieinsum!(ein"ij,jk->ik", (a, b), c)
	@test c ≈ a*b
end
