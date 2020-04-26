function contract4224!(::CC{DEVICE,1}, A, B, C, D) where {DEVICE}
    ein"abcd,(be,(hc,efgh))->afgd"(A, B, C, D)
end

function contract223a(cc::CC{DEVICE,1}, A, B, C) where {DEVICE}
    ein"(ab,bc),cde->ade"(A, B, C)
end

function contract222a(cc::CC{DEVICE,1}, A, B, C) where {DEVICE}
    ein"ab,bc,cd->ad"(A, B, C)
end

function contract222b(cc::CC{DEVICE,1}, A, B, C) where {DEVICE}
    ein"ab,da,cd->cb"(A, B, C)
end

function contract223b(cc::CC{DEVICE,1}, A, B, C) where {DEVICE}
    ein"(ab,ae),cde->cdb"(A, B, C)
end

function contract4223(cc::CC{DEVICE,1}, A, B, C, D) where {DEVICE}
    ein"abcd,(be,(hc,efh))->afd"(A, B, C, D)
end

function contract3223(cc::CC{DEVICE,1}, A, B, C, D) where {DEVICE}
    ein"abc,(ad,(fb,def))->ec"(A, B, C, D)
end

function contract3222(cc::CC{:CPU}, A, B, C, D) where {DEVICE}
    ein"abc,(ad,(fb,df))->c"(A, B, C, D)
end

function contract3222(cc::CC{:GPU,1}, A, B, C, D) where {DEVICE}
    ein"abc,(ad,(fb,df))->c"(A, B, C, D)
end

function contract3222(cc::CC{:GPU,N}, A, B, C, D) where N
    ein"abc,(ad,(fb,df))->c"(A, Array(B), Array(C), Array(D))
end

array_upload(::CC{:GPU}, A, inds) = CuArray(A[inds...])
array_upload(::CC{:CPU}, A, inds) = Array(A[inds...])
function array_download!(A, resi, inds)
    A[inds...] = resi
end

function contract4224!(cc::CC{DEVICE,N}, A, B, C, D) where {N,DEVICE}
    L = size(A, 1)
    M = size(A, 4)
    _con!(cc, code)
end

function _con!(cc, code, out, A, Ts...; dim::Int)
    L = size(A, dim)
    chunk_size = ceil(Int, L ÷ N)
    for i = 1:N
        if chunk_size*(i-1) < min(L,chunk_size*i)
            slice = Any[Colon() for i=1:ndims(A)]
            slice[dim] = chunk_size*(i-1)+1:min(L,chunk_size*i)
            Ai = array_upload(cc, A, slice)
            # using ein"abcd,be,hc,efgh->afgd" is very slow, seems like a bug.
            resi = code(Ai, Ts...)
            array_download!(out, resi, slice)
        end
    end
    out
end
