using TropicalTensors, BenchmarkTools, CUDA
using LinearAlgebra, TropicalNumbers
using BenchmarkTools
using CUDA
CUDA.allowscalar(false)

function LinearAlgebra.mul!(C::Matrix{T}, A::Matrix{T}, B::Matrix{T}) where T<:TropicalTypes
    @assert size(A, 2) == size(B, 1) && size(A, 1) == size(C, 1) && size(B, 2) == size(C, 2)
    for j=1:size(B, 2)
        for i=1:size(A, 1)
            for k=1:size(A, 2)
                @inbounds C[i,j] += A[i, k] * B[k, j]
            end
        end
    end
    C
end

a = randn(1000, 1000)
b = randn(1000, 1000)
cua = CuArray(a);
cub = CuArray(b);
@benchmark a * b
@benchmark cua * cub

a = randn(CountingTropical{Float64}, 1000, 1000)
b = randn(CountingTropical{Float64}, 1000, 1000)
cua = CuArray(a);
cub = CuArray(b);
@benchmark a * b
@benchmark cua * cub

A = zeros(CountingTropical{Float64}, 10, 32, 21);
B = zeros(CountingTropical{Float64}, 32, 11, 5, 23, 41, 10);
#tensorcontract((1,2,3), A, (2,4,5,6,7,1), B, (7,4,3,5,6))
@benchmark tensorcontract((1,2,3), A, (2,4,5,6,7,1), B, (7,4,3,5,6))

cuA = CuArray(A);
cuB = CuArray(B);
@benchmark CUDA.@sync tensorcontract((1,2,3), cuA, (2,4,5,6,7,1), cuB, (7,4,3,5,6))