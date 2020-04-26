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

function _con!(cc::CC{DEVICE,1}, code, out, A, Ts...; dim::Int) where {DEVICE}
    copyto!(out, code(A, Ts...))
end

_con(cc::CC{DEVICE,1}, code, A, Ts...; dim::Int) where {DEVICE} = code(A, Ts...)

function _con!(cc::CC{DEVICE,N}, code, out, A, Ts...; dim::Int) where {DEVICE,N}
    L = size(A, dim)
    ix = get_leftmost_ixs(code)
    odim = indexin(ix[dim], [OMEinsum.getiy(code.eins)...])[]
    chunk_size = ceil(Int, L/N)
    for i = 1:N
        if chunk_size*(i-1) < min(L,chunk_size*i)
            range = chunk_size*(i-1)+1:min(L,chunk_size*i)
            slice = Any[Colon() for i=1:ndims(A)]
            slice[dim] = range
            Ai = array_upload(cc, A, slice)
            # using ein"abcd,be,hc,efgh->afgd" is very slow, seems like a bug.
            resi = code(Ai, Ts...)
            slice = Any[Colon() for i=1:ndims(out)]
            slice[odim] = range
            array_download!(out, resi, slice)
        end
    end
    out
end

function _con(cc::CC{DEVICE,N}, code, A, Ts...; dim::Int) where {DEVICE,N}
    iy = OMEinsum.getiy(code.eins)
    isize = get_size(code, (A, Ts...), iy)
    out = similar(A, isize)
    _con!(cc::CC{DEVICE,N}, code, out, A, Ts...; dim=dim)
    out
end

macro split(cc, dim, expr::Expr)
    :($_con($cc, $(expr.args[1]), $(expr.args[2:end]...); dim=$dim)) |> esc
end

macro split!(cc, dim, expr::Expr)
    :($_con!($cc, $(expr.args[1]), $(expr.args[2]), $(expr.args[2:end]...); dim=$dim)) |> esc
end

function get_leftmost_ixs(code::OMEinsum.NestedEinsumStable)
    l = code.args[1]
    if l isa Int
        OMEinsum.getixs(code.eins)[1]
    else
        get_leftmost_ixs(l)
    end
end

function get_size(code::OMEinsum.NestedEinsumStable, xs, iy)
    ixl, ixr = OMEinsum.getixs(code.eins)
    lsize = get_size(code.args[1], xs, ixl)
    rsize = get_size(code.args[2], xs, ixr)
    d = OMEinsum.IndexSize((ixl..., ixr...), (lsize..., rsize...))
    map(i->d[i], iy)
end

function get_size(code::Int, xs, iy)
    xs[code] |> size
end
