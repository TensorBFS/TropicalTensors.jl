# Solving the Potts model on a SquareLattice

export potts_bondtensor, δ, solve_potts

"""
    potts_bondtensor(T, Val(q), J; digits=5)

`T` is the matrix element type,
`J` is the coupling strength,
`q` is the degree of freedom of each spin,
`digits` will make the floating point numbers to have better numeric properties, it is important to the success of degeneracy counting!
"""
function potts_bondtensor(::Type{T}, ::Val{q}, J; digits=5) where {T, q}
    angles = cos.(2π .* ((1:q) ./ q))
    res = zeros(T, q, q)
    for i=1:q
        for j=1:q
            res[i,j] = T(round(J*angles[mod1(abs(j-i), q)]; digits=5))
        end
    end
    res
end

function δ(::Type{T}, ::Val{q}, n::Int) where {T, q}
    res = zeros(T, fill(q, n)...)
    for i=1:q
        res[fill(i, n)...] = one(T)
    end
    res
end

function sequential_tree(n::Int)
    n == 1 && return 1
    TropicalTensors.ContractionTree(sequential_tree(n-1), n)
end

function solve_potts(::Type{T}, ::Val{q}, lt::SquareLattice, J::Dict; usecuda=false, digits=5) where {T, q}
    # generate tensors
    tensors = Matrix{Any}(undef,lt.Nx, lt.Ny)

    for i=1:lt.Nx
        for j=1:lt.Ny
            v = δ(T, Val(q), 4)
            labels = [(i,j-1)=>(i,j), (i-1,j)=>(i,j), (i,j)=>(i,j+1), (i,j)=>(i+1,j)]
            mask = [j==1, i==1, j==lt.Ny, i==lt.Nx]
            if i!=lt.Nx
                b = potts_bondtensor(T, Val(q), J[(i,j)=>(i+1,j)]; digits=digits)
                v = tensorcontract((1,2,3,4), v, (4,5), b, (1,2,3,5))
            end
            if j!=lt.Ny
                a = potts_bondtensor(T, Val(q), J[(i,j)=>(i,j+1)]; digits=digits)
                v = tensorcontract((1,2,3,4), v, (3,5), a, (1,2,5,4))
            end
            w = dropdims(sum(v, dims=findall(mask)); dims=(findall(mask)...,))
            tensors[i,j] = LabeledTensor(w, labels[(!).(mask)])
        end
    end
    tn = TensorNetwork(tensors)
    if usecuda
        tn = TropicalTensors.togpu(tn)
    end

    # construct contraction tree
    tree = sequential_tree(length(tensors))

    # contract
    Array(TropicalTensors.contract(tn, tree).array)[]
end

function potts_randpm_J(lt::SquareLattice)
    d = Dict{Pair{Tuple{Int,Int},Tuple{Int,Int}},Int}()
    for i=1:lt.Nx
        for j=1:lt.Ny
            if i!=lt.Nx
                d[(i,j)=>(i+1,j)] = rand([1,-1])
            end
            if j!=lt.Ny
                d[(i,j)=>(i,j+1)] = rand([1,-1])
            end
        end
    end
    return d
end

