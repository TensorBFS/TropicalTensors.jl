using EliminateGraphs
using TropicalTensors
using OMEinsum
using CuArrays
CuArrays.allowscalar(false)

struct NoSplit end
# split first dimension
struct Split{N} end

function get_TroI(::Type{Array})
    I2 = zeros(2,2)
    I2 .= TROPICAL_ZERO
    I2[1,1] = TROPICAL_ONE
    I2[2,2] = TROPICAL_ONE

    I3 = zeros(2,2,2)
    I3 .= TROPICAL_ZERO
    I3[1,1,1] = TROPICAL_ONE
    I3[2,2,2] = TROPICAL_ONE

    I4 = zeros(2,2,2,2)
    I4 .= TROPICAL_ZERO
    I4[1,1,1,1] = TROPICAL_ONE
    I4[2,2,2,2] = TROPICAL_ONE

    return Tropical.(I2),Tropical.(I3),Tropical.(I4)
end

_get_J(::Val{:ferro}) = -1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()

function get_TroB(::Type{Array}, v::Val)
    Jij = _get_J(v)
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

get_TroB(::Type{CuArray}, v::Val) = CuArray(get_TroB(Array, v))
get_TroI(::Type{CuArray}) = CuArray.(get_TroI(Array))

function gen_mine_2d(L::Int, Jtype::Val, array_type::Type, contract_method)
    I2,I3,I4 = get_TroI(array_type)
    """The first row"""
    ball = I2
    for j in 2:L-1
        ball = ein"(ab,bc),cde->ade"(ball,get_TroB(array_type, Jtype),I3)
        ball = reshape( ball,:,size(ball,3) )
    end
    ball = ein"ab,bc,cd->ad"(ball,get_TroB(array_type, Jtype),I2)

    """ From the second row to the second row from bottom """
    for i in 2:L-1
        ball = reshape(ball,2,:)
        ball = ein"(ab,ae),cde->cdb"(ball,get_TroB(array_type, Jtype),I3)
        ball = reshape(ball,2,2,2,:)
        for j in 2:L-1
            ball = contract4224!(contract_method, ball,get_TroB(array_type, Jtype),get_TroB(array_type, Jtype),I4)
            ball = reshape(ball, size(ball,1)*size(ball,2),2,2,:)
        end
        ball = ein"abcd,(be,(hc,efh))->afd"(ball,get_TroB(array_type, Jtype),get_TroB(array_type, Jtype),I3)
        ball = reshape(ball,2,:)
    end

    """ the last row """
    ball = ein"ab,da,cd->cb"(ball,get_TroB(array_type, Jtype),I2)
    ball = reshape(ball,2,2,:)
    for j in 2:L-1
        ball = ein"abc,(ad,(fb,def))->ec"(ball,get_TroB(array_type, Jtype),get_TroB(array_type, Jtype),I3)
        ball = reshape(ball,2,2,:)
    end
    ball = ein"abc,(ad,(fb,df))->c"(ball,get_TroB(array_type, Jtype),get_TroB(array_type, Jtype),I2)
    return ball
end

function contract4224!(::NoSplit, A::AT, B, C, D) where AT
    ein"abcd,(be,(hc,efgh))->afgd"(A, B, C, D)
end

array_upload(::Type{<:CuArray}, A) = CuArray(A)
array_upload(::Type{<:Array}, A) = Array(A)

function contract4224!(::Split{N}, A::AT, B, C, D::DT) where {N,AT,DT}
    L = size(A, 1)
    M = size(A, 4)
    chunk_size = ceil(Int, (L >= M ? L : M) ÷ N)
    for i = 1:N
        if L >= M
            if chunk_size*(i-1) < min(L,chunk_size*i)
                Ai = array_upload(DT, view(A, chunk_size*(i-1)+1:min(L,chunk_size*i),:,:,:))
                # using ein"abcd,be,hc,efgh->afgd" is very slow, seems like a bug.
                resi = ein"abcd,(be, (hc,efgh))->afgd"(Ai, B, C, D)
                A[chunk_size*(i-1)+1:min(L,chunk_size*i),:,:,:] = resi
            end
        else
            if chunk_size*(i-1) < min(M,chunk_size*i)
                # TODO: fix zero dimension error!
                Ai = array_upload(DT, view(A, :,:,:,chunk_size*(i-1)+1:min(M,chunk_size*i)))
                resi = ein"abcd,(be, (hc,efgh))->afgd"(Ai, B, C, D)
                A[:,:,:,chunk_size*(i-1)+1:min(M,chunk_size*i)] = resi
            end
        end
    end
    A
end

L = 16
e1 = gen_mine_2d(L, Val(:ferro), CuArray, Split{4}())
e2 = gen_mine_2d(L, Val(:ferro), Array, Split{4}())
e3 = gen_mine_2d(L, Val(:ferro), CuArray, NoSplit())
e4 = gen_mine_2d(L, Val(:ferro), Array, NoSplit())
using Test
@test Array(e1) == e2 == Array(e3) == e4

using BenchmarkTools
@benchmark gen_mine_2d($L, $(Val(:ferro)), $Array)

using CuArrays
e = gen_mine_2d(L, Val(:ferro), CuArray)
@benchmark (CuArrays.@sync gen_mine_2d($L, $(Val(:ferro)), $CuArray))
