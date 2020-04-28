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

    return Tropical.(Float32.(I2)),Tropical.(Float32.(I3)),Tropical.(Float32.(I4))
end

function get_TroI64(::Type{Array})
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

    return Tropical.(Float64.(I2)),Tropical.(Float64.(I3)),Tropical.(Float64.(I4))
end
_get_J(::Val{:ferro}) = 1.0
_get_J(::Val{:anti}) = -1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()
_get_J(::Val{:pm}) = rand() > 0.5 ? 1.0 : -1.0

function get_TroB(::Type{Array}, v::Val)
    Jij = Float32(_get_J(v))
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

function get_TroB64(::Type{Array}, v::Val)
    Jij = Float64(_get_J(v))
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

#get_TroB(::Type{CuArray}, v::Val) = CuArray(get_TroB(Array, v))
#get_TroI(::Type{CuArray}) = CuArray.(get_TroI(Array))

function gen_2d(L::Int,array_type::Type, Jtype::Val)
    # [L rows, L columns] each of which contains L-1 edges
    return [[[get_TroB64(array_type,Jtype) for j in 1:L-1] for i in 1:L ] for k in 1:2]
end

function get_T(i::Int,j::Int, L::Int, array_type::Type, config::Array{Array{Int,1},1})
    I2,I3,I4 = get_TroI64(array_type)
    T = i == 1 || i == L ? (T = j==1 || j==L ? I2 : I3) : (T = j==1 || j==L ? I3 : I4)
    if config[i][j] == 1
        T[end] = Tropical(-Inf)
    elseif config[i][j] == 2
        T[1] = Tropical(-Inf)
    end

    return T
end

function noninf_ratio(ball)
    #println("# ",sum(ball .== Tropical(-Inf))," -Inf, totally ", length(ball))
    return
end

function mine_2d(L::Int, Jtype::Val, array_type::Type, Js, config::Array{Array{Int,1},1} )
    """
    return minimum energy of the 2D model.
    """
    I2,I3,I4 = get_TroI(array_type)
    rows,cols = Js
    """The first row"""
    i=1
    j=1
    ball = get_T(i,j,L,array_type,config)
    #println("# ",1.0-sum(ball .== Tropical(-Inf)) / length(ball))
    noninf_ratio(ball)
    for j in 2:L-1
        #ball = ein"ab,bc,cde->ade"(ball,rows[i][j-1],I3)
        ball = ein"ab,bc,cde->ade"(ball,rows[i][j-1],get_T(i,j,L,array_type,config))
        ball = reshape( ball,:,size(ball,3) )
    end
    j=L
    #println("# ",1.0-sum(ball .== Tropical(-Inf)) / length(ball))
    noninf_ratio(ball)
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
        #println("# ",1.0-sum(ball .== Tropical(-Inf)) / length(ball))
        noninf_ratio(ball)
    end

    """ the last row """
    i = L
    j = 1
    #ball = ein"ab,da,cd->cb"(ball,cols[j][i-1],I2)
    ball = ein"ab,da,cd->cb"(ball,cols[j][i-1],get_T(i,j,L,array_type,config))
    ball = reshape(ball,2,2,:)
    #println("# ",1.0-sum(ball .== Tropical(-Inf)) / length(ball))
    noninf_ratio(ball)
    for j in 2:L-1
        #ball = ein"abc,ad,fb,def->ec"(ball,rows[i][j-1],cols[j][i-1],I3)
        ball = ein"abc,ad,fb,def->ec"(ball,rows[i][j-1],cols[j][i-1],get_T(i,j,L,array_type,config))
        ball = reshape(ball,2,2,:)
    end
    j = L
    #ball = ein"abc,ad,fb,df->c"(ball,rows[j][i-1],cols[j][i-1],I2)
    ball = ein"abc,ad,fb,df->c"(ball,rows[j][i-1],cols[j][i-1],get_T(i,j,L,array_type,config))
    println(typeof(ball))
    return ball[1]
end

function gen_2d_2(L::Int,array_type::Type, Jtype::Val)
    # [L rows, L columns] each of which contains L edges, first of which is 1.
    return [[ pushfirst!([get_TroB(array_type,Jtype) for j in 1:L-1],Tropical.(zeros(Float32,1,1))) for i in 1:L ] for k in 1:2]
end

function Js_(Js)
    #return [[[ pushfirst!(k, Tropical.(zeros(1,1))) for k in j] for j in i  ] for i in Js]
    return [[ pushfirst!(j, Tropical.(zeros(Float32,1,1)))  for j in i  ] for i in Js]
end

function get_T_2(i::Int,j::Int, L::Int, array_type::Type, config::Array{Array{Int,1},1},Js)
    @assert 1<=i<=L && 1<=j<=L
    rows,cols = Js
    I2,I3,I4 = get_TroI(array_type)
    T = i == 1 || i == L ? (T = j==1 || j==L ? I2 : I3) : (T = j==1 || j==L ? I3 : I4)
    if config[i][j] == 1
        T[end] = Tropical(-Inf32)
    elseif config[i][j] == 2
        T[1] = Tropical(-Inf32)
    end

    if i==1
        if j==1
            T = reshape(T,1,2,2,1)
        elseif j<L
            T = reshape(T,2,2,2,1)
        else
            T = reshape(T,2,2,1,1)
        end
    elseif i<L
        if j==1
            T = reshape(T,1,2,2,2)
        elseif j<L
            T = reshape(T,2,2,2,2)
        else
            T = reshape(T,2,2,1,2)
        end
    else
        if j==1
            T = reshape(T,1,1,2,2)
        elseif j<L
            T = reshape(T,2,1,2,2)
        else
            T = reshape(T,2,1,1,2)
        end
    end
    #println("(",i,")",j," ",size(rows[i][j]),size(cols[j][i]),size(T))
    return ein"ab,bcde,fe->acdf"(rows[i][j],T,cols[j][i])
end

function mine_2d_2(L::Int, Jtype::Val, array_type::Type, Js, config::Array{Array{Int,1},1} )
    """
    return minimum energy of the 2D model.
    """
    ball = Tropical.(zeros(Float32,1,1))
    for i in 1:L
        ball = reshape(ball,1,:)
        for j in 1:L
            T = get_T_2(i,j,L,array_type,config,Js)
            ball = ein"abcd,bfgc->afgd"( reshape(ball, size(ball,1),size(T,1),size(T,4),:) , T)
            ball = reshape(ball, size(ball,1)*size(ball,2),:)
        end
    end
    println(typeof(ball))
    return ball[1]
end

function mine_2d_gpu(L::Int, Jtype::Val, Js, config::Array{Array{Int,1},1},split::Int )
    """
    return minimum energy of the 2D model.
    """
    ball = Tropical.(zeros(Float32,1,1))
    for i in 1:L
        ball = reshape(ball,1,:)
        for j in 1:L
            T = get_T_2(i,j,L,array_type,config,Js)
            ball = ein"abcd,bfgc->afgd"( reshape(ball, size(ball,1),size(T,1),size(T,4),:) , T)
            ball = reshape(ball, size(ball,1)*size(ball,2),:)
        end
    end
    println(typeof(ball))
    return ball[1]
end

function ground_state(L::Int, Jtype::Val, array_type::Type,Js)
    config = [ [0 for j in 1:L] for i in 1:L ]
    println("energy of the found state is ",mine_2d(L, Jtype, Array, Js, config))
    for i in 1:L
        for j in 1:L
            config[i][j] = 1
            ep = mine_2d(L, Jtype, Array,Js, config)
            config[i][j] = 2
            em = mine_2d(L, Jtype, Array, Js,config)
            #println("i= ",i," j=",j," ep=",ep," em=",em)
            if ep > em
                config[i][j] = 1
            else
                config[i][j] = 2
            end
        end
    end
    println("Ground state: ",config)
    println("Its energy is ",mine_2d(L, Jtype, Array, Js, config))
end


L = 16
#Jtype = Val(:ferro)
#Jtype = Val(:pm)
#Jtype = Val(:randn)
#Js = gen_2d(L,Array,Jtype)
#config = [ [0 for j in 1:L] for i in 1:L ]
#ground_state(L, Jtype, Array,Js)
println("L= ",L," Jtype",Jtype)
#@time e = mine_2d(L, Jtype, Array,Js, config)

#Js2 = Js_(Js)
Js2 = gen_2d_2(L,Array,Jtype)
@time e = mine_2d_2(L, Jtype, Array,Js2, config)
#@time e = mine_2d(L, Val(:anti), Array, config)
#@time e = mine_2d(L, Val(:randn), Array, config)
#println("mine=",e)
