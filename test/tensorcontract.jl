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

@testset "tensor contract" begin
    A = zeros(Float64, 10, 32, 21);
    B = zeros(Float64, 32, 11, 5, 2, 41, 10);
    tA = LabeledTensor(A, [1,2,3])
    tB = LabeledTensor(B, [2,4,5,6,7,1])

    tOut = tA * tB
    @test tOut.array == ein"abc,bdefga->cdefg"(A, B)
    @test tOut.labels == [3,4,5,6,7]

    tnet1 = TensorNetwork([tA])
    @test tnet1 isa TensorNetwork
    tnet = TensorNetwork([tA, tB])

    tOut2, contracted_labels = contract_label!(tnet, 1)
    @test tnet.tensors[] ≈ tOut
    @test tnet.tensors[] ≈ tOut2
    @test contracted_labels == [1, 2]
end
