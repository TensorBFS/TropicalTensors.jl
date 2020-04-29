using TropicalTensors
using OMEinsum
using CuArrays
using Printf
CuArrays.allowscalar(false)

function get_TroI(dtype)
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
    if dtype == Float32
        return Tropical.(Float32.(I2)),Tropical.(Float32.(I3)),Tropical.(Float32.(I4))
    elseif dtype == Float16
        return Tropical.(Float16.(I2)),Tropical.(Float16.(I3)),Tropical.(Float16.(I4))
    end
end


function get_TroI(array_type,dtype)
    if array_type == CuArray
        return CuArray.(get_TroI(dtype))
    elseif array_type == Array
        return get_TroI(dtype)
    end
end


_get_J(::Val{:ferro}) = 0.1
_get_J(::Val{:anti}) = -1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()
_get_J(::Val{:pm}) = rand() > 0.5 ? 1.0 : -1.0

function get_TroB(dtype, v::Val)
    if dtype == Float32
        Jij = Float32(_get_J(v))
    elseif dtype == Float16
        Jij = Float16(_get_J(v))
    end
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

function get_TroB(array_type,dtype, v::Val)
    if array_type == CuArray
        return CuArray(get_TroB(dtype,v))
    elseif array_type == Array
        return get_TroB(dtype,v)
    end
end

function gen_2d(L::Int,atype,dtype, Jtype::Val)
    # [L rows, L columns] each of which contains L edges, first of which is 1.
    if atype == Array
        return [[ pushfirst!([get_TroB(atype,dtype,Jtype) for j in 1:L-1],Tropical.(zeros(dtype,1,1))) for i in 1:L ] for k in 1:2]
    elseif atype == CuArray
        return [[ pushfirst!([get_TroB(CuArray,dtype,Jtype) for j in 1:L-1],(CuArray(Tropical.(zeros(dtype,1,1))))) for i in 1:L ] for k in 1:2]
    end
end

function get_T(i::Int,j::Int, L::Int, atype::Type, Js,dtype::Type)
    @assert 1<=i<=L && 1<=j<=L
    rows,cols = Js
    I2,I3,I4 = get_TroI(atype,dtype)
    T = i == 1 || i == L ? (T = j==1 || j==L ? I2 : I3) : (T = j==1 || j==L ? I3 : I4)

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

function mine_2d(L::Int, Jtype::Val, atype::Type, Js, dtype )
    """
    return minimum energy of the 2D model.
    """
    ball = Tropical.(zeros(dtype,1,1))
    if atype == CuArray
        ball = CuArray(ball)
    end
    for i in 1:L
        ball = reshape(ball,1,:)
        for j in 1:L
            t0 = time()
            print("\r(",i,",",j,") / (",L,",",L,") \t")
            T = get_T(i,j,L,atype,Js,dtype)
            n_chunks = 4
            if i==1 || j==1 || j==L || i==L
                ball = ein"abcd,bfgc->afgd"( reshape(ball, size(ball,1),size(T,1),size(T,4),:) , T)
            elseif j <= L/2
                ball = reshape(ball, size(ball,1),size(T,1),size(T,n_chunks),:,n_chunks)
                for chunk in 1:n_chunks
                    ball[:,:,:,:,chunk] = ein"abcd,bfgc->afgd"( ball[:,:,:,:,chunk] , T)
                end
            else
                ball = reshape(ball,n_chunks,Int(size(ball,1)/n_chunks),size(T,1),size(T,4),:)
                for chunk in 1:n_chunks
                    ball[chunk,:,:,:,:] = ein"abcd,bfgc->afgd"( ball[chunk,:,:,:,:] , T)
                end
                ball = reshape(ball,size(ball,1)*size(ball,2),size(T,2),:)
            end
            ball = reshape(ball, size(ball,1)*size(ball,2),:)
            @printf("%.2g Sec.                  ",time()-t0)
        end
    end
    println(" Done")
    return ball
end




L=32
Jtype = Val(:ferro)
#dtype = Float32
dtype = Float16
#atype = Array
atype = CuArray
Js = gen_2d(L,atype, dtype, Jtype)
@time mine=mine_2d(L, Jtype, atype, Js, dtype )
println("mine = ",mine)

