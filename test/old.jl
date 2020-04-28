using EliminateGraphs
using TropicalTensors
using OMEinsum
using CuArrays
CuArrays.allowscalar(false)

function get_TroI(::Type{Array})
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

function get_TroB(::Type{Array}, v::Val)
    Jij = _get_J(v)
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

get_TroB(::Type{CuArray}, v::Val) = CuArray(get_TroB(Array, v))
get_TroI(::Type{CuArray}) = CuArray.(get_TroI(Array))

function gen_mine_2d(L::Int, Jtype::Val, array_type::Type)
    I2,I3,I4 = get_TroI(array_type)
    """The first row"""
    ball = I2
    for j in 2:L-1
        ball = ein"ab,bc,cde->ade"(ball,get_TroB(array_type, Jtype),I3)
        ball = reshape( ball,:,size(ball,3) )
    end
    ball = ein"ab,bc,cd->ad"(ball,get_TroB(array_type, Jtype),I2)

    """ From the second row to the second row from bottom """
    for i in 2:L-1
        ball = reshape(ball,2,:)
        ball = ein"ab,ae,cde->cdb"(ball,get_TroB(array_type, Jtype),I3)
        ball = reshape(ball,2,2,2,:)
        for j in 2:L-1
            ball = ein"abcd,be,hc,efgh->afgd"(ball,get_TroB(array_type, Jtype),get_TroB(array_type, Jtype),I4)
            ball = reshape(ball, size(ball,1)*size(ball,2),2,2,:)
        end
        ball = ein"abcd,be,hc,efh->afd"(ball,get_TroB(array_type, Jtype),get_TroB(array_type, Jtype),I3)
        ball = reshape(ball,2,:)
    end

    """ the last row """
    ball = ein"ab,da,cd->cb"(ball,get_TroB(array_type, Jtype),I2)
    ball = reshape(ball,2,2,:)
    for j in 2:L-1
        ball = ein"abc,ad,fb,def->ec"(ball,get_TroB(array_type, Jtype),get_TroB(array_type, Jtype),I3)
        ball = reshape(ball,2,2,:)
    end
    ball = ein"abc,ad,fb,df->c"(ball,get_TroB(array_type, Jtype),get_TroB(array_type, Jtype),I2)
    return ball
end

L = 16
e = gen_mine_2d(L, Val(:ferro), Array)
println("mine=",e)

using BenchmarkTools
@benchmark gen_mine_2d($L, $(Val(:ferro)), $Array)

using CuArrays
e = gen_mine_2d(L, Val(:ferro), CuArray)
@benchmark (CuArrays.@sync gen_mine_2d($L, $(Val(:ferro)), $CuArray))
