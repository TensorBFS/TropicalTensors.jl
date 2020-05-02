using TropicalTensors
using Base.Threads
using OMEinsum
using Base.Cartesian

for T in [:Float16, :Float32, :Float64]
    elsize = if T == :Float16
        16
    elseif T == :Float32
        8
    elseif T == :Float64
        4
    end
    B = 32 ÷ elsize
    @eval function tropical_mm!(c::StridedMatrix{Tropical{$T}}, a::StridedMatrix{Tropical{$T}}, b::StridedMatrix{Tropical{$T}})
        @assert (size(a, 1) ÷ $B) * $B == size(a, 1)
        @assert size(a, 2) == size(b, 1)
        @inbounds @threads for k=1:size(a, 2)
            offset_a = $elsize*size(a, 1)*(k-1)
            for j=1:size(b,2)
                offset_c = (j-1)*size(c, 1)*$elsize
                bb = b[k,j].n
                for i=1:size(a, 1)÷$B
                    offset_ai = offset_a + 32 * i
                    offset_ci = offset_c + 32 * i
                    ai = vload(SVec{$B,$T}, pointer(a) + offset_ai)
                    ai = vadd(ai, bb)
                    c[i,j] = Tropical(max(c[i,j].n, maximum(ai)))
                    #vstore!(pointer(a) + offset, ai)
                    #@nexprs $B il->if ci < ai[il]
                        #c[offset_ci,j].n = Tropical(ab1)
                    #end
                end
            end
        end
        return c
    end
end

function tropical_mm2!(c::StridedMatrix{Tropical{T}}, a::StridedMatrix{Tropical{T}}, b::StridedMatrix{Tropical{T}}) where T
    @assert size(a, 2) == size(b, 1)
    @inbounds @threads for k=1:size(a, 2)
        for j=1:size(b,2)
            bb = b[k,j].n
            for i=1:size(a, 1)
                ab = a[i,k].n + bb
                ci = c[i,j].n
                if ab > ci
                    c[i,j] = Tropical(ab)
                end
            end
        end
    end
    return c
end

a = Tropical.(randn(512, 512))
b = Tropical.(randn(512, 512))
c = Tropical.(zeros(512, 512))

using BenchmarkTools
@benchmark ein"ab,bc->ac"($a, $b)
@benchmark (c = zero($a); tropical_mm!(c, $a, $b))
@benchmark (c = zero($a); tropical_mm2!(c, $a, $b))
tropical_mm!(zero(a), a, b) ≈ a * b
tropical_mm2!(zero(a), a, b) ≈ a * b
