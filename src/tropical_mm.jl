using Base.Threads
using LinearAlgebra

function LinearAlgebra.mul!(c, a::AbstractMatrix{<:Tropical}, b::AbstractMatrix{<:Tropical})
    @assert size(a, 2) == size(b, 1)
    for k=1:size(a, 2)
        @threads for j=1:size(b,2)
            @inbounds @simd for i=1:size(a, 1)
                ab = a[i,k].n + b[k,j].n
                if c[i,j].n < ab
                    c[i,j] = Tropical(ab)
                end
            end
        end
    end
    return c
end
