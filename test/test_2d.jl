using EliminateGraphs
using TropicalTensors
using OMEinsum
using CuArrays
CuArrays.allowscalar(false)

# Computing Configuration
# N is the number of chunks
struct CC{DEVICE,N} end

function get_TroI(::CC{:CPU})
    I2 = zeros(2,2)
    I2 .= -Inf
    I2[1,1] = 0.0
    I2[2,2] = 0.0

    I3 = zeros(2,2,2)
    I3 .= -Inf
    I3[1,1,1] = 0.0
    I3[2,2,2] = 0.0

    I4 = zeros(2,2,2,2)
    I4 .= -Inf
    I4[1,1,1,1] = 0.0
    I4[2,2,2,2] = 0.0

    return Tropical.(I2),Tropical.(I3),Tropical.(I4)
end

_get_J(::Val{:ferro}) = 1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()

function get_TroB(::CC{:CPU}, v::Val)
    Jij = _get_J(v)
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

get_TroB(::CC{:GPU}, v::Val) = CuArray(get_TroB(Array, v))
get_TroI(::CC{:GPU}) = CuArray.(get_TroI(Array))

function gen_mine_2d(L::Int, Jtype::Val, ::CC)
    I2,I3,I4 = get_TroI(cc)
    """The first row"""
    ball = I2
    for j in 2:L-1
        ball = ein"(ab,bc),cde->ade"(ball,get_TroB(cc, Jtype),I3)
        ball = reshape( ball,:,size(ball,3) )
    end
    ball = ein"ab,bc,cd->ad"(ball,get_TroB(cc, Jtype),I2)

    """ From the second row to the second row from bottom """
    for i in 2:L-1
        ball = reshape(ball,2,:)
        ball = ein"(ab,ae),cde->cdb"(ball,get_TroB(cc, Jtype),I3)
        ball = reshape(ball,2,2,2,:)
        for j in 2:L-1
            ball = contract4224!(cc, ball,get_TroB(cc, Jtype),get_TroB(cc, Jtype),I4)
            ball = reshape(ball, size(ball,1)*size(ball,2),2,2,:)
        end
        ball = ein"abcd,(be,(hc,efh))->afd"(ball,get_TroB(cc, Jtype),get_TroB(cc, Jtype),I3)
        ball = reshape(ball,2,:)
    end

    """ the last row """
    ball = ein"ab,da,cd->cb"(ball,get_TroB(cc, Jtype),I2)
    ball = reshape(ball,2,2,:)
    for j in 2:L-1
        ball = ein"abc,(ad,(fb,def))->ec"(ball,get_TroB(cc, Jtype),get_TroB(cc, Jtype),I3)
        ball = reshape(ball,2,2,:)
    end
    ball = ein"abc,(ad,(fb,df))->c"(ball,get_TroB(cc, Jtype),get_TroB(cc, Jtype),I2)
    return ball
end

function contract4224!(::NoSplit, A::AT, B, C, D) where AT
    ein"abcd,(be,(hc,efgh))->afgd"(A, B, C, D)
end

array_upload(::CC{:GPU}, A) = CuArray(A)
array_upload(::CC{:CPU}, A) = Array(A)

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
e1 = gen_mine_2d(L, Val(:ferro), CC{:GPU,4}())
e2 = gen_mine_2d(L, Val(:ferro), CC{:CPU,4}())
e3 = gen_mine_2d(L, Val(:ferro), CC{:GPU,1}())
e4 = gen_mine_2d(L, Val(:ferro), CC{:CPU,1}())
using Test
@test Array(e1) == e2 == Array(e3) == e4

using BenchmarkTools
@benchmark gen_mine_2d($L, $(Val(:ferro)), CC{:CPU,1}())

using CuArrays
e = gen_mine_2d(L, Val(:ferro), CC{:GPU}())
@benchmark (CuArrays.@sync gen_mine_2d($L, $(Val(:ferro)), CC{:GPU}()))
