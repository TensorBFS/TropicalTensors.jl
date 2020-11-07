using Test
using OMEinsum
using TropicalTensors

@testset "contract" begin
    for (iAs, iBs, iOut) in [
        [(1,2,3), (2,3,4), (1,4)],
        [(1,3,2), (2,3,4), (4,1)],
        [(), (2,4), (4,2)],
        [(2,4), (), (4,2)],
        [(2,4), (1,), (4,2,1)],
        [(2,4), (2,4), ()],
        [(2,4), (4,), (2,)],
        ]
        A = asarray(randn(rand(4:12, length(iAs))...))
        B = asarray(randn([(iBs[i] in iAs) ? size(A, indexin(iBs[i], collect(iAs))[]) : rand(4:12) for i=1:length(iBs)]...))
        out = tensorcontract(iAs, A, iBs, B, iOut)
        eincode = EinCode{(iAs, iBs), iOut}()
        eout = eincode(A, B)
        @test eout ≈ out
    end
end