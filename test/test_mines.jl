"""
Find ground state energy and entropy
The couplings are integers, so that there are degeneracies.
"""
function logsumexp(vec::Array{Float64,1})
    maxval = maximum(vec)
    return log(sum(exp.(vec .- maxval)))+maxval
end

function tro_mm(ball::Tuple{Array{Int64,2},Array{Float64,2}},B::Array{Int64,2},my_half_inf::Int64)
    """
    Tropical Matrix Multiplication with log counting
    A,S = ball where A stores the energy information and S stores the entropy information
    """
    A,S = ball
    n,m = size(A)
    n2,m2 = size(S)
    k,l = size(B)
    @assert m==k && n==n2 && m==m2
    C = repeat( reshape(A,n,1,m),1,l,1 ) .+ repeat( reshape(B',1,l,k),n,1,1 )
    D = reshape(maximum(C,dims=3),n,l)
    S_expand = repeat( reshape(S,n,1,m),1,l,1 ) # (n,l,m)
    S = zeros(n,l)
    for i in 1:n
        for j in 1:l
            idx = (C[i,j,:] .== D[i,j])
            S[i,j] = logsumexp(S_expand[i,j,idx])
        end
    end
    D[D .< -my_half_inf] .= -my_half_inf*2
    return Tuple((D,S))
end

function get_I(::Type{Array},my_half_inf::Int)
    I2 = Int.(zeros(2,2))
    I2 .= -my_half_inf*2
    I2[1,1] = 0
    I2[2,2] = 0

    I3 = Int.(zeros(2,2,2))
    I3 .= -my_half_inf*2
    I3[1,1,1] = 0
    I3[2,2,2] = 0

    I4 = Int.(zeros(2,2,2,2))
    I4 .= -my_half_inf*2
    I4[1,1,1,1] = 0
    I4[2,2,2,2] = 0
    return I2,I3,I4
end

_get_J(::Val{:ferro}) = Int(1)
_get_J(::Val{:anti}) = -Int(1)
_get_J(::Val{:pm}) = rand() > 0.5 ? Int(1) : -Int(1)


function get_B(::Type{Array}, v::Val)
    Jij = _get_J(v)
    return [Jij -Jij; Jij -Jij]
end

function gen_2d(L::Int,array_type::Type, Jtype::Val)
    # [L rows, L columns] each of which contains L-1 edges
    return [[[get_B(array_type,Jtype) for j in 1:L-1] for i in 1:L ] for k in 1:2]
end

function get_T(i::Int,j::Int, L::Int, array_type::Type, config::Array{Array{Int,1},1},my_half_inf::Int)
    I2,I3,I4 = get_I(array_type,my_half_inf)
    T = i == 1 || i == L ? (T = j==1 || j==L ? I2 : I3) : (T = j==1 || j==L ? I3 : I4)
    if config[i][j] == 1
        T[end] = -my_half_inf*2
    elseif config[i][j] == 2
        T[1] = -my_half_inf*2
    end
    return T
end

function mines_2d(L::Int, Jtype::Val, array_type::Type, Js, config::Array{Array{Int64,1},1} )
    """
    return minimum energy of the 2D model.
    The negative infinity is maintained by hand using "Int(2*my_half_inf)", because Inf is Float...
        my_half_inf = upperbound(-E),
        my_inf = 2*my_half_inf, to make sure that my_inf + upperbound(-E) will be still recognizable.
    """
    my_half_inf = Int(sum(sum([ sum([ sum([(abs.(reshape(k,1,:))) for k in j]) for j in i]) for i in Js]))/4 + 1)
    I2,I3,I4 = get_I(array_type,my_half_inf)
    rows,cols = Js
    """The first row"""
    i = 1
    j = 1
    T = reshape(get_T(i,j,L,array_type,config,my_half_inf),:,2)
    ball = Tuple((T, zeros(size(T))))
    for j in 2:L
        ball = tro_mm(ball,rows[i][j-1],my_half_inf)
        myI = get_T(i,j,L,array_type,config,my_half_inf)
        ball = tro_mm(ball,reshape(myI,2,:),my_half_inf)
        ball = Tuple(( reshape(A,:,2) for A in ball ))
    end
    """ From the second row to the last row """
    for i in 2:L
        j=1
        ball = tro_mm( Tuple(( permutedims(reshape(A,2,:),[2,1]) for A in ball )),cols[j][i-1],my_half_inf )
        ball = tro_mm( Tuple(ball), reshape(get_T(i,j,L,array_type,config,my_half_inf),2,:),my_half_inf)
        for j in 2:L
            ball = tro_mm( Tuple((reshape(A,:,2) for A in ball)), rows[i][j-1],my_half_inf)
            ball = ( permutedims(reshape(A,2,:,2),[2,3,1]) for A in ball )
            ball = tro_mm( Tuple((reshape(A,:,2) for A in ball )),cols[j][i-1],my_half_inf)
            ball = ( reshape(A,:,4) for A in ball )
            ball = tro_mm(Tuple(ball),reshape(get_T(i,j,L,array_type,config,my_half_inf),4,:),my_half_inf)
        end
    end
    return ball
end

L = 10
#Jtype = Val(:ferro)
#Jtype = Val(:anti)
Jtype = Val(:pm)
Js = gen_2d(L,Array,Jtype)
config = [ [0 for j in 1:L] for i in 1:L ]
E,S = mines_2d(L, Jtype, Array, Js, config )
println("E= ",-E[1]," \t Entropy= ",S[1])
