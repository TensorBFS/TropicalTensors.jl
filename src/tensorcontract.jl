export tensorcontract

function _align_eltypes(xs::AbstractArray...)
    T = promote_type(eltype.(xs)...)
    return map(x->eltype(x)==T ? x : T.(x), xs)
end

function _align_eltypes(xs::AbstractArray{T}...) where T
    xs
end

# batched routines
@inline _indexpos(iAs, i)::Int = findfirst(==(i), iAs)
_isunique(x) = length(unique(x)) == length(x)

# can be used in either static or dynamic invoke
@noinline function analyse_dim(iAs, iBs, iOuts)
    # check indices
    @assert _isunique(iAs) "indices in A matrix is not unique $iAs"
    @assert _isunique(iBs) "indices in B matrix is not unique $iBs"
    @assert _isunique(iOuts) "indices in C matrix is not unique $iOuts"
    allinds = [iAs..., iBs..., iOuts...]
    @assert all(i -> count(==(i), allinds) == 2, allinds) "indices does not appear in pairs! got $iAs, $iBs and $iOuts"

    mdims = setdiff(iAs, iBs)
    ndims = setdiff(iBs, iAs)
    kdims = iAs ∩ iBs
    iABs = mdims ∪ ndims
    lA = mdims ∪ kdims
    lB = kdims ∪ ndims
    lOut = mdims ∪ ndims
    pA = indexin(lA, iAs)
    pB = indexin(lB, iBs)
    pOut = indexin(iOuts, lOut)
    pA, pB, pOut, length(kdims)
end

@noinline function analyse_size(pA, sA, pB, sB, nc)
    nA = length(sA)
    nB = length(sB)
    sA1 = Int[sA[pA[i]] for i=1:nA-nc]
    sB2 = Int[sB[pB[i]] for i=nc+1:nB]
    K = mapreduce(i->sB[pB[i]], *, 1:nc, init=1)
    return prod(sA1), K, prod(sB2), [sA1..., sB2...]
end

function tensorcontract(iAs, A::AbstractArray, iBs, B::AbstractArray, iOuts)
    A, B = _align_eltypes(A, B)
    pA, pB, pOut, nc = analyse_dim([iAs...], [iBs...], [iOuts...])
    M, K, N, sOut = analyse_size(pA, size(A), pB, size(B), nc)

    Apr = reshape(_conditioned_permutedims(A, pA), M, K)
    Bpr = reshape(_conditioned_permutedims(B, pB), K, N)
    AB = Apr * Bpr
    AB = _conditioned_permutedims(reshape(AB, sOut...), (pOut...,))
end

function _conditioned_permutedims(A::AbstractArray{T,N}, perm) where {T,N}
    if any(i-> (@inbounds perm[i]!=i), 1:N)
        return permutedims(A, perm)
    else
        return A
    end
end