export hypercubicI, copyvertex, copytensor, copygate
export Gcp, Greset
export MaskedSquareLattice, rand_maskedsquare

struct MaskedSquareLattice <: Viznet.AbstractSquareLattice
    mask::Matrix{Bool}
end

Base.size(lt::MaskedSquareLattice) = size(lt.mask)
Base.size(lt::MaskedSquareLattice, i::Int) = size(lt.mask, i)
Viznet.unit(lt::MaskedSquareLattice) = 1/(max(size(lt)...))
Viznet.vertices(lt::MaskedSquareLattice) = [i for i in 1:length(lt.mask) if lt.mask[i]]
function rand_maskedsquare(Nx::Int, Ny::Int, ρ::Real; seed=2)
    Random.seed!(seed)
    MaskedSquareLattice(rand(Nx, Ny) .< ρ)
end
function Base.getindex(lt::MaskedSquareLattice, i::Real, j::Real)
    step = unit(lt)
    (i-0.5)*step, (j-0.5)*step
end

function Viznet.isconnected(sq::MaskedSquareLattice, i::Int, j::Int)
    u = unit(sq)
    d = distance(sq[i], sq[j])
    sq.mask[i] && sq.mask[j] && d > 0.999u && d < (sqrt(2)+0.001)*u
end

Viznet.bonds(lt::MaskedSquareLattice) = sgbonds(lt)

function solve(::Type{TT}, sg::Spinglass{LT,T}; usecuda=false) where {TT, LT<:MaskedSquareLattice,T}
    lt = sg.lattice
    Lx, Ly = size(lt)
    nbit = Ly + 2
    reg = _init_reg(TT, lt, Val(usecuda))
    Js = copy(sg.Js)
    hs = copy(sg.hs)
    LI = LinearIndices(lt)
    _c(a, b) = isconnected(lt, LI[a...], LI[b...])
    for i=1:Lx
        println("Layer $i/$Lx")
        for j=1:Ly-1
            _c((i,j), (i,j+1)) && (reg |> put(nbit, (j,j+1)=>Gvb(bondtensor(TT, sg, Js |> popfirst!))))
        end
        for j=1:Ly
            if lt.mask[i,j]
                reg |> put(nbit, j=>Gh(vertextensor(TT, sg, hs |> popfirst!)))
            else
                reg |> put(nbit, j=>Gcut(TT))
            end
        end
        (i!=Lx) && for j=1:Ly
            ancpre = nbit-(j-1)%2
            ancthis = nbit-(j)%2
            # store the information in qubit `j` to ancilla `nbit-j%2`
            j!=Ly && _c((i,j), (i+1,j+1)) && (reg |> put(nbit, (j, ancthis)=>Gcp(TT)))
            # interact with j-1 th qubit (a)
            j!=1 && _c((i+1,j-1), (i,j)) && (reg |> put(nbit, (j-1,j)=>Gvb(bondtensor(TT, sg, Js |> popfirst!))))
            # onsite term (b)
            _c((i,j), (i+1,j)) && (reg |> put(nbit, j=>Ghb(bondtensor(TT, sg, Js |> popfirst!))))
            if j!=1 && _c((i,j-1), (i+1,j))
                # interact with cached j-1 th qubit (c)
                jj = Js[1]
                reg |> put(nbit, (ancpre,j)=>Gvb(bondtensor(TT, sg, Js |> popfirst!)))
                # erease the information in previous ancilla
                reg |> put(nbit, ancpre=>Gcut(TT))
            end
        end
    end
    sum(state(reg))
end

function sgvertices(lt::MaskedSquareLattice)
    [v for v in sgvertices(SquareLattice(size(lt)...)) if lt.mask[v]]
end

function sgbonds(lt::MaskedSquareLattice)
    edges = Tuple{Int, Int}[]
    LI = LinearIndices(size(lt))
    function trypush!(i, j)
        if lt.mask[i...] && lt.mask[j...]
            push!(edges, (LI[i...], LI[j...]))
        end
    end
    Lx, Ly = size(lt)
    for i=1:Lx
        for j=1:Ly-1
            trypush!([i,j], [i,j+1])
        end
        (i!=Lx) && for j=1:Ly
            j!=1 && trypush!([i+1,j-1], [i,j])
            trypush!([i,j], [i+1,j])
            j!=1 && trypush!([i,j-1],[i+1,j])
        end
    end
    edges
end

regsize(lt::MaskedSquareLattice) = size(lt,2)+2

function assign_Js_hs(lt::MaskedSquareLattice, grad_Js::AbstractVector{T}, grad_hs) where T
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

function _init_reg(::Type{T}, lt::MaskedSquareLattice, ::Val{:false}) where T
    nbit = size(lt, 2) + 2
    state = ones(T, 1<<nbit)
    ArrayReg(state)
end
