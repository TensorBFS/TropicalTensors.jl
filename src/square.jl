export sgbonds

function solve(::Type{TT}, sg::AbstractSpinglass{LT}; usecuda=false) where {LT<:SquareLattice,TT}
    # Yao gates
    lt = sg.lattice
    Lx, Ly = lt.Nx, lt.Ny
    reg = _init_reg(TT, lt.Ny, Val(usecuda))
    k = 0

    println("Layer 1/$Lx")
    for j=1:Ly
        reg |> put(Ly, j=>Gh(vertextensor(TT, sg, j)))
    end
    for j=1:Ly-1
        k += 1
        reg |> put(Ly, (j,j+1)=>Gvb(bondtensor(TT, sg, k)))
    end
    for i=2:Lx
        println("Layer $i/$Lx")
        for j=1:Ly
            k += 1
            reg |> put(Ly, j=>Ghb(bondtensor(TT, sg, k)))
        end
        for j=1:Ly
            reg |> put(Ly, j=>Gh(vertextensor(TT, sg, (i-1)*Ly + j)))
        end
        for j=1:Ly-1
            k += 1
            reg |> put(Ly, (j,j+1)=>Gvb(bondtensor(TT, sg, k)))
        end
    end
    sum(state(reg))
end

function solve(lattice::Viznet.AbstractLattice, Js::Vector{T}, hs::Vector{T}; usecuda=false) where {T}
    sg = Spinglass(lattice, Js, hs)
    solve(sg; usecuda=usecuda)
end

# to index `h`.
function sgvertices(lt::SquareLattice)
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
cachesize_A(lt::SquareLattice) = lt.Ny
cachesize_B(lt::SquareLattice) = lt.Nx-1
cachesize_largemem(lt::SquareLattice) = (lt.Nx-1) * lt.Ny

function _init_reg(::Type{T}, lt::SquareLattice, usecuda) where T<:TropicalTypes
    _init_reg(T, lt.Ny, usecuda)
end
