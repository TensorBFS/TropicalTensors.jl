using EliminateGraphs
using TropicalTensors
using OMEinsum
#using CuArrays
#CuArrays.allowscalar(false)

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
_get_J(::Val{:anti}) = 1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()

function get_TroB(::Type{Array}, v::Val)
    Jij = _get_J(v)
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

#get_TroB(::Type{CuArray}, v::Val) = CuArray(get_TroB(Array, v))
#get_TroI(::Type{CuArray}) = CuArray.(get_TroI(Array))

function gen_2d(L::Int,array_type::Type, Jtype::Val)
    # [L rows, L columns] each of which contains L-1 edges
    return [[[get_TroB(array_type,Jtype) for j in 1:L-1] for i in 1:L ] for k in 1:2]
end

function get_T(i::Int,j::Int, L::Int, array_type::Type, config::Array{Int64,2})
    I2,I3,I4 = get_TroI(array_type)
    if i == 1 || i == L
        if j == 1 || j == L
            T = I2
        else
            T = I3
        end
    else
        if j == 1 || j == L
            T = I3
        else
            T = I4
        end
    end
    if j == 1 || j == L
        if i == 1 || i == L
            T = I2
        else
            T = I3
        end
    else
        if i == 1 || i == L
            T = I3
        else
            T = I4
        end
    end

    if config[i,j] == 1
        T[end] = TROPICAL_ZERO
    elseif config[i,j] == 2
        T[1] = TROPICAL_ZERO
    end

    return T
end

function mine_2d(LA::Int, Jtype::Val, array_type::Type, Js, config::Array{Int, 2} )
    """
    return minimum energy of the 2D model.
    """
    I2,I3,I4 = get_TroI(array_type)
    rows,cols = Js
    """The first row"""
    i=1
    j=1
    ball = get_T(i,j,L,array_type,config)
    for j in 2:L-1
        #ball = ein"ab,bc,cde->ade"(ball,rows[i][j-1],I3)
        ball = ein"ab,bc,cde->ade"(ball,rows[i][j-1],get_T(i,j,L,array_type,config))
        ball = reshape( ball,:,size(ball,3) )
    end
    j=L
    #ball = ein"ab,bc,cd->ad"(ball,rows[i][j-1],I2)
    ball = ein"ab,bc,cd->ad"(ball,rows[i][j-1],get_T(i,j,L,array_type,config))

    """ From the second row to the second row from bottom """
    for i in 2:L-1
        ball = reshape(ball,2,:)
        j=1
        #ball = ein"ab,ae,cde->cdb"(ball,cols[j][i-1],I3)
        ball = ein"ab,ae,cde->cdb"(ball,cols[j][i-1],get_T(i,j,L,array_type,config))
        ball = reshape(ball,2,2,2,:)
        for j in 2:L-1
            #ball = ein"abcd,be,hc,efgh->afgd"(ball,rows[i][j-1],cols[j][i-1],I4)
            ball = ein"abcd,be,hc,efgh->afgd"(ball,rows[i][j-1],cols[j][i-1],get_T(i,j,L,array_type,config))
            ball = reshape(ball, size(ball,1)*size(ball,2),2,2,:)
        end
        j = L
        #ball = ein"abcd,be,hc,efh->afd"(ball,rows[i][j-1],cols[j][i-1],I3)
        ball = ein"abcd,be,hc,efh->afd"(ball,rows[i][j-1],cols[j][i-1],get_T(i,j,L,array_type,config))
        ball = reshape(ball,2,:)
    end

    """ the last row """
    i = L
    j = 1
    #ball = ein"ab,da,cd->cb"(ball,cols[j][i-1],I2)
    ball = ein"ab,da,cd->cb"(ball,cols[j][i-1],get_T(i,j,L,array_type,config))
    ball = reshape(ball,2,2,:)
    for j in 2:L-1
        #ball = ein"abc,ad,fb,def->ec"(ball,rows[i][j-1],cols[j][i-1],I3)
        ball = ein"abc,ad,fb,def->ec"(ball,rows[i][j-1],cols[j][i-1],get_T(i,j,L,array_type,config))
        ball = reshape(ball,2,2,:)
    end
    j = L
    #ball = ein"abc,ad,fb,df->c"(ball,rows[j][i-1],cols[j][i-1],I2)
    ball = ein"abc,ad,fb,df->c"(ball,rows[j][i-1],cols[j][i-1],get_T(i,j,L,array_type,config))
    return ball[1]
end

function ground_state(L::Int, Jtype::Val, array_type::Type,Js)
    config = zeros(Int, L, L)
    println("energy of the found state is ",mine_2d(L, Jtype, Array, Js, config))
    for i in 1:L
        for j in 1:L
            config[i,j] = 1
            ep = mine_2d(L, Jtype, Array,Js, config)
            config[i,j] = 2
            em = mine_2d(L, Jtype, Array, Js,config)
            #println("i= ",i," j=",j," ep=",ep," em=",em)
            if ep < em
                config[i,j] = 1
            else
                config[i,j] = 2
            end
        end
    end
    println("Ground state: ",config)
    println("Its energy is ",mine_2d(L, Jtype, Array, Js, config))
end


L = 10
#Jtype = Val(:ferro)
#Jtype = Val(:anti)
Jtype = Val(:randn)
Js = gen_2d(L,Array,Jtype)
ground_state(L, Jtype, Array,Js)

#config = [ [0 for j in 1:L] for i in 1:L ]
#@time e = mine_2d(L, Val(:ferro), Array, config)
#@time e = mine_2d(L, Val(:anti), Array, config)
#@time e = mine_2d(L, Val(:randn), Array, config)
#println("mine=",e)
