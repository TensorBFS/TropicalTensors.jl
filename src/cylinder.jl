export hypercubicI, copyvertex, copytensor, copygate
export Cylinder, rand_maskedsquare

"""periodic in y direction."""
struct Cylinder <: Viznet.AbstractSquareLattice
    Nx::Int
    Ny::Int
end

Base.size(lt::Cylinder) = (lt.Nx, lt.Ny)
Base.size(lt::Cylinder, i::Int) = size(lt, i)
Viznet.unit(lt::Cylinder) = 1/(max(size(lt)...))
Viznet.vertices(lt::Cylinder) = 1:lt.Nx * lt.Ny
function Base.getindex(lt::Cylinder, i::Real, j::Real)
    step = unit(lt)
    (i-0.5)*step, (j-0.5)*step
end

function solve(::Type{TT}, sg::AbstractSpinglass{LT}; usecuda=false) where {LT<:Cylinder,TT}
    # Yao gates
    lt = sg.lattice
    Lx, Ly = lt.Nx, lt.Ny
    reg = _init_reg(TT, lt.Ny, Val(usecuda))
    k = 0

    for i=1:Lx
        println("Layer $i/$Lx")
        i!=1 && for j=1:Ly
            k += 1
            reg |> put(Ly, j=>Ghb(bondtensor(TT, sg, k)))
        end
        for j=1:Ly
            reg |> put(Ly, j=>Gh(vertextensor(TT, sg, (i-1)*Ly + j)))
        end
        for j=1:Ly
            k += 1
            reg |> put(Ly, (j,mod1(j+1, Ly))=>Gvb(bondtensor(TT, sg, k)))
        end
    end
    sum(state(reg))
end

Viznet.bonds(lt::Cylinder) = sgbonds(lt)

function sgvertices(lt::Cylinder)
    sgvertices(SquareLattice(size(lt)...))
end

function sgbonds(lt::Cylinder)
    edges = Tuple{Int, Int}[]
    LI = LinearIndices(size(lt))
    function trypush!(i, j)
        push!(edges, (LI[i...], LI[j...]))
    end
    Lx, Ly = size(lt)
    for i=1:Lx
        for j=1:Ly
            trypush!([i,j], [i,mod1(j+1,Ly)])
        end
        (i!=Lx) && for j=1:Ly
            trypush!([i,j], [i+1,j])
        end
    end
    edges
end

function assign_Js_hs(lt::Cylinder, grad_Js::AbstractVector{T}, grad_hs) where T
    grid = zeros(size(lt.mask))
    vorder = sgvertices(lt)
    for (i, hg) in enumerate(grad_hs)
        grid[vorder[i]] = hg
    end
    for ((i, j), Jg) in zip(sgbonds(lt), grad_Js)
        assign_one!(grid, i, j, Jg)
    end
    return grid
end

regsize(lt::Cylinder) = size(lt,2)
function _init_reg(::Type{T}, lt::Cylinder, ::Val{:false}) where T
    nbit = regsize(lt)
    state = ones(T, 1<<nbit)
    ArrayReg(state)
end
