using TropicalTensors

function potts3_bondtensor(::Type{T}, J, ::Val{q}) where {T, q}
    angles = cos.(2π .* ((1:q) ./ q))
    res = zeros(T, q, q)
    for i=1:q
        for j=1:q
            res[i,j] = T(J*angles[mod1(abs(j-i), q)])
        end
    end
    res
end

function δ3(::Type{T}, n::Int) where T
    res = zeros(T, fill(3, n)...)
    for i=1:3
        res[fill(i, n)...] = one(T)
    end
    res
end

function sequential_tree(n::Int)
    n == 1 && return 1
    TropicalTensors.ContractionTree(sequential_tree(n-1), n)
end

function solve_potts3(::Type{T}, lt::SquareLattice, J::Dict; usecuda=false) where T
    # generate tensors
    tensors = Matrix{Any}(undef,lt.Nx, lt.Ny)

    for i=1:lt.Nx
        for j=1:lt.Ny
            v = δ3(T, 4)
            labels = [(i,j-1)=>(i,j), (i-1,j)=>(i,j), (i,j)=>(i,j+1), (i,j)=>(i+1,j)]
            mask = [j==1, i==1, j==lt.Ny, i==lt.Nx]
            if i!=lt.Nx
                b = potts3_bondtensor(T, J[(i,j)=>(i+1,j)], Val(3))
                v = tensorcontract((1,2,3,4), v, (4,5), b, (1,2,3,5))
            end
            if j!=lt.Ny
                a = potts3_bondtensor(T, J[(i,j)=>(i,j+1)], Val(3))
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
    TropicalTensors.contract_tree(tn, tree).array[].n
end

function build_J(lt)
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

lt = SquareLattice(9, 9)
res = solve_potts3(Tropical{Float64}, lt, build_J(lt); usecuda=false)