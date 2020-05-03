using SIMDPirates
using TropicalTensors
using Base.Threads

for T in [:Float16, :Float32, :Float64]
    elsize = if T == :Float16
        2
    elseif T == :Float32
        4
    elseif T == :Float64
        8
    end
    B = 32 ÷ elsize
    @eval function tropical_mm_simd!(c::StridedMatrix{Tropical{$T}}, a::StridedMatrix{Tropical{$T}}, b::StridedMatrix{Tropical{$T}})
        @assert (size(a, 1) ÷ $B) * $B == size(a, 1)
        @assert size(a, 2) == size(b, 1)
        @inbounds for k=1:size(a, 2)
            offset_a = $elsize*size(a, 1)*(k-1)
            @threads for j=1:size(b,2)
                offset_c = (j-1)*size(c, 1)*$elsize
                bb = b[k,j].n
                for i=1:size(a, 1)÷$B
                    offset_ai = offset_a + 32 * (i-1)
                    offset_ci = offset_c + 32 * (i-1)
                    ai = vload(SVec{$B,$T}, pointer(a) + offset_ai)
                    ci = vload(SVec{$B,$T}, pointer(c) + offset_ci)
                    ai = vadd(ai, bb)
                    vstore!(Ptr{$T}(pointer(c) + offset_ci), max(ai, ci))
                end
            end
        end
        return c
    end
end

using Test
@testset "mm" begin
    a = Tropical.(randn(512, 512))
    b = Tropical.(randn(512, 512))
    @test tropical_mm_simd!(zero(a), a, b) ≈ a * b
end
