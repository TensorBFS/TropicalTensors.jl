using OMEinsum
using TropicalNumbers

function get_TroI(::Type{T}, N) where T
    I2 = zeros(T, fill(2, N)...)
    I2[1] = one(T)
    I2[end] = one(T)
    return I2
end

function get_TroB(Jij::T) where T
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

maxdim(A) = findmax(size(A))[2]
function gen_mine_2d(L::Int, Js::AbstractVector{T}) where T
    I2, I3, I4 = get_TroI.(Tropical{T}, (2, 3, 4))
    Js = copy(Js)
    # The first row
    ball = I2
    for j in 2:L-1
        ball = ein"(ab,bc),cde->ade"(ball, get_TroB(popfirst!(Js)), I3)
        ball = reshape(ball,:,size(ball,3))
    end
    ball = ein"(ab,bc),cd->ad"(ball,get_TroB(popfirst!(Js)), I2)

    # From the second row to the second row from bottom
    for i in 2:L-1
        ball = reshape(ball,2,:)
        ball = ein"(ab,ae),cde->cdb"(ball, get_TroB(popfirst!(Js)), I3)
        ball = reshape(ball,2,2,2,:)
        for j in 2:L-1
            ball = ein"abcd,(be,(hc,efgh))->afgd"(ball, get_TroB(popfirst!(Js)), get_TroB(popfirst!(Js)), I4)
            ball = reshape(ball, size(ball,1)*size(ball,2),2,2,:)
        end
        ball = ein"abcd,(be,(hc,efh))->afd"(ball, get_TroB(popfirst!(Js)),get_TroB(popfirst!(Js)), I3)
        ball = reshape(ball,2,:)
    end

    # the last row
    ball = ein"(ab,da),cd->cb"(ball,get_TroB(popfirst!(Js)), I2)
    ball = reshape(ball, 2, 2, :)
    for j in 2:L-1
        ball = ein"abc,(ad,(fb,def))->ec"(ball, get_TroB(popfirst!(Js)), get_TroB(popfirst!(Js)), I3)
        ball = reshape(ball, 2, 2, :)
    end
    ball = ein"abc,(ad,(fb,df))->c"(ball, get_TroB(popfirst!(Js)),get_TroB(popfirst!(Js)), I2)
    return ball
end

L = 10
gen_mine_2d(L, ones(L * (L-1) * 2))
