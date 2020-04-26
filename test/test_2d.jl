using EliminateGraphs
using TropicalTensors
using OMEinsum
using CuArrays
CuArrays.allowscalar(false)

Base.isapprox(x::Tropical, y::Tropical; kwargs...) = Base.isapprox(x.n, y.n; kwargs...)
Base.isapprox(x::Array{<:Tropical}, y::Array{<:Tropical}; kwargs...) = all(Base.isapprox.(x, y; kwargs...))
Base.rtoldefault(x::Type{Tropical{Float64}}, ::Type{Tropical{Float64}}, ::Int64)  where T = 1e-7
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

get_TroB(::CC{:GPU}, v::Val) = CuArray(get_TroB(CC{:CPU,1}(), v))
get_TroI(::CC{:GPU}) = CuArray.(get_TroI(CC{:CPU,1}()))

include("contractions.jl")

maxdim(A) = findmax(size(A))[2]
function gen_mine_2d(L::Int, Jtype::Val, cc::CC)
    I2,I3,I4 = get_TroI(cc)
    CC1 = CC{:CPU, 1}()
    # The first row
    ball = I2
    for j in 2:L-1
        ball = @split cc 1 ein"(ab,bc),cde->ade"(ball,get_TroB(cc, Jtype),I3)
        ball = reshape(ball,:,size(ball,3))
    end
    ball = @split cc 1 ein"(ab,bc),cd->ad"(ball,get_TroB(cc, Jtype),I2)

    # From the second row to the second row from bottom
    for i in 2:L-1
        ball = reshape(ball,2,:)
        ball = @split cc 2 ein"(ab,ae),cde->cdb"(ball,get_TroB(cc, Jtype),I3)
        ball = reshape(ball,2,2,2,:)
        for j in 2:L-1
            ball = @split! cc maxdim(ball) ein"abcd,(be,(hc,efgh))->afgd"(ball,get_TroB(cc, Jtype),get_TroB(cc, Jtype),I4)
            ball = reshape(ball, size(ball,1)*size(ball,2),2,2,:)
        end
        ball = @split cc maxdim(ball) ein"abcd,(be,(hc,efh))->afd"(ball,get_TroB(cc, Jtype),get_TroB(cc, Jtype),I3)
        ball = reshape(ball,2,:)
    end

    # the last row
    ball = @split cc 2 ein"(ab,da),cd->cb"(ball,get_TroB(cc, Jtype),I2)
    ball = reshape(ball,2,2,:)
    for j in 2:L-1
        ball = @split cc 3 ein"abc,(ad,(fb,def))->ec"(ball,get_TroB(cc, Jtype),get_TroB(cc, Jtype),I3)
        ball = reshape(ball,2,2,:)
    end
    ball = @split cc 3 ein"abc,(ad,(fb,df))->c"(ball,get_TroB(cc, Jtype),get_TroB(cc, Jtype),I2)
    return ball
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
e = gen_mine_2d(L, Val(:ferro), CC{:GPU,1}())
@benchmark (CuArrays.@sync gen_mine_2d($L, $(Val(:ferro)), CC{:GPU,1}()))
