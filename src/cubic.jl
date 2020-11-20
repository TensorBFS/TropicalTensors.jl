export CubicLattice

struct CubicLattice <: Viznet.AbstractLattice
    Nx::Int
    Ny::Int
    Nz::Int
end

Base.length(lt::CubicLattice) = lt.Nx * lt.Ny * lt.Nz

function getbonds(lt::CubicLattice)
    out = Pair{Tuple{Int,Int,Int}, Tuple{Int,Int,Int}}[]
    for k=1:lt.Nz, j=1:lt.Ny, i=1:lt.Nx
        i!=lt.Nx && push!(out, (i,j,k)=>(i+1,j,k))
        j!=lt.Ny && push!(out, (i,j,k)=>(i,j+1,k))
        k!=lt.Nz && push!(out, (i,j,k)=>(i,j,k+1))
    end
    out
end

function getvertices(lt::CubicLattice)
    out = Tuple{Int,Int,Int}[]
    for k=1:lt.Nz, j=1:lt.Ny, i=1:lt.Nx
        push!(out, (i,j,k))
    end
    out
end

function rand_Js(::Type{T}, ::Ferro, lt::CubicLattice) where T
    Dict(b=>one(T) for b in getbonds(lt))
end
function rand_Js(::Type{T}, ::Randpm, lt::CubicLattice) where T
    Dict(b=>T(rand([1,-1])) for b in getbonds(lt))
end
function rand_hs(::Type{T}, ::Zero, lt::CubicLattice) where T
    Dict(b=>zero(T) for b in getvertices(lt))
end

function solve(::Type{TT}, sg::AbstractSpinglass{LT}; usecuda=false) where {LT<:CubicLattice,TT}
    # Yao gates
    lt = sg.lattice
    Lx, Ly, Lz = lt.Nx, lt.Ny, lt.Nz
    LI = LinearIndices((Lx, Ly))
    n = Lx * Ly
    reg = _init_reg(TT, n, Val(usecuda))

    for k=1:Lz
        println("Layer $k/$Lz")
        k!=1 && for j=1:Ly, i=1:Lx
            reg |> put(n, LI[i,j]=>Ghb(bondtensor(TT, sg, (i,j,k-1)=>(i,j,k))))
        end
        for j=1:Ly, i=1:Lx
            reg |> put(n, LI[i,j]=>Gh(vertextensor(TT, sg, (i,j,k))))
        end
        for j=1:Ly, i=1:Lx
            j!=Ly && (reg |> put(n, (LI[i,j],LI[i,j+1])=>Gvb(bondtensor(TT, sg, (i,j,k)=>(i,j+1,k)))))
            i!=Lx && (reg |> put(n, (LI[i,j],LI[i+1,j])=>Gvb(bondtensor(TT, sg, (i,j,k)=>(i+1,j,k)))))
        end
    end
    sum(state(reg))
end

regsize(lt::CubicLattice) = lt.Ny * lt.Nx

function _init_reg(::Type{T}, lt::CubicLattice, usecuda) where T<:TropicalTypes
    _init_reg(T, regsize(lt), usecuda)
end
