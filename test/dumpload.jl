using TropicalTensors
using Test

@testset "dump load" begin
    jt = Randn()
    ht = Randn()
    rsg = rand_spinglass(Float64, SquareLattice(3,3); jt=jt, ht=ht, dumpto=".")
    rsg2 = rand_spinglass(Float64, SquareLattice(3,3))
    load_params(".", rsg2, jt, ht, 2)
    @test rsg.Js ≈ rsg2.Js
    @test rsg.hs ≈ rsg2.hs
end
