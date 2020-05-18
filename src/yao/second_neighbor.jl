export hypercubicI, copyvertex, copytensor, copygate
export Gcp, Greset
export SquareLattice2nd

struct SquareLattice2nd <: Viznet.AbstractSquareLattice
    Nx::Int
    Ny::Int
end

Base.size(lt::SquareLattice2nd) = (lt.Nx, lt.Ny)
Viznet.vertices(lt::SquareLattice2nd) = 1:lt.Nx*lt.Ny

function sgbonds(lt::SquareLattice2nd)
    edges = Tuple{Int, Int}[]
    LI = LinearIndices(size(lt))
    Lx, Ly = size(lt)
    for i=1:Lx
        for j=1:Ly-1
            push!(edges, (LI[i,j], LI[i,j+1]))
        end
        (i!=Lx) && for j=1:Ly
            j!=1 && push!(edges, (LI[i+1,j-1], LI[i,j]))
            push!(edges, (LI[i,j], LI[i+1,j]))
            j!=1 && push!(edges, (LI[i,j-1],LI[i+1,j]))
        end
    end
    edges
end

function sgvertexorder(lt::SquareLattice2nd)
    sgvertexorder(SquareLattice(lt.Nx, lt.Ny))
end

"""
    hypercubicI([T], ndim::Int, D::Int)

Get a `ndim`-dimensional identity hypercubic, with bound dimension `D`.
"""
function hypercubicI(::Type{T}, ndim::Int, D::Int) where T
    res = zeros(T, fill(D, ndim)...)
    @inbounds for i=1:D
        res[fill(i,ndim)...] = one(T)
    end
    return res
end

hypercubicI(ndim::Int, D::Int) = hypercubicI(Float64, ndim, D)
function copyvertex(::Type{T}) where T
    res = zeros(T, 2, 2, 2)
    res[1,1,1] = one(T)
    res[2,2,1] = one(T)
    res
end

function copytensor(::Type{T}) where T
    PermMatrix([1,4,3,2], [one(T), zero(T), zero(T), one(T)])
end

"""
    copygate(T)

copy state of qubit 2 -> 1.
"""
function Gcp(::Type{T}) where T
    matblock(copytensor(Tropical{T}))
end

function Greset(::Type{T}) where T
    TT = Tropical{T}
    matblock([one(TT) one(TT); zero(TT) zero(TT)])
end

function solve(sg::Spinglass{LT,T}; usecuda=false) where {LT<:SquareLattice2nd,T}
    lt = sg.lattice
    Lx, Ly = size(lt)
    nbit = Ly + 2
    reg = _init_reg(T, nbit, Val(usecuda))
    Js = copy(sg.Js)
    hs = copy(sg.hs)
    for i=1:Lx
        println("Layer $i/$Lx")
        for j=1:Ly-1
            reg |> put(nbit, (j,j+1)=>G4(T, Js |> popfirst!))
        end
        for j=1:Ly
            reg |> put(nbit, j=>Gh(T, hs |> popfirst!))
        end
        (i!=Lx) && for j=1:Ly
            # store the information in qubit `j` to ancilla `nbit-j%2`
            j!=Ly && (reg |> put(nbit, (j,nbit-j%2)=>Gcp(T)))
            # interact with j-1 th qubit (a)
            j!=1 && (reg |> put(nbit, (j-1,j)=>G4(T, Js |> popfirst!)))
            # onsite term (b)
            reg |> put(nbit, j=>G2(T, Js |> popfirst!))
            # interact with cached j-1 th qubit (c)
            j!=1 && (reg |> put(nbit, (nbit-(j-1)%2,j)=>G4(T, Js |> popfirst!)))
            # erease the information in previous ancilla `nbit-(j-1)%2`
            j!=1 && (reg |> put(nbit, nbit-(j-1)%2=>Greset(T)))
        end
    end
    sum(state(reg))
end

# to index `h`.
function sgvertexorder(lt::SquareLattice)
    v = Int[]
    LI = LinearIndices(lt)
    for i=1:lt.Nx
        for j=1:lt.Ny
            push!(v, LI[i,j])
        end
    end
    return v
end

# to index `J`.
function sgbonds(lt::SquareLattice)
    edges  =Tuple{Int,Int}[]
    append!(edges, v_bonds(lt, 1))
    for i=2:lt.Nx
        append!(edges, h_bonds(lt, i-1))
        append!(edges, v_bonds(lt, i))
    end
    edges
end

function v_bonds(lt::SquareLattice, i::Int)
    LI = LinearIndices(size(lt))
    map(j -> (LI[i,j], LI[i,j+1]), 1:lt.Ny-1)
end

function h_bonds(lt::SquareLattice, i::Int)
    LI = LinearIndices(size(lt))
    map(j -> (LI[i,j], LI[i+1,j]), 1:lt.Ny)
end

regsize(lt::SquareLattice) = lt.Ny
