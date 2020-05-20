function red_reg(::Type{T}, Ly::Int, Js, hs; usecuda=false) where T
    reg = _init_reg(T, Ly*4, Val(usecuda))
    k = 0
    for i=1:Ly*4
        hs[i] != 0 && (reg |> put(4*Ly, i=>Gh(T, hs[i])))
    end
    for i=1:Ly-1
        for j = 1:4
            reg |> put(4*Ly, (4*(i-1)+j,4*i+j)=>G4(T, Js[k+1]))
            k += 1
        end
    end
    for i=1:Ly
        reg |> put(4*Ly, (4i-3,4i-2,4i-1,4i)=>G16(T, Js[k+1:k+16]))
        k += 16
    end
    @assert k==length(Js)
    return reg
end

function solve(sg::Spinglass{LT,T}; usecuda::Bool) where {LT<:ChimeraLattice,T}
    _, Lx, Ly = size(sg.lattice)
    Js = sg.Js
    hs = sg.hs
    println("Running Chimera, Lx = $Lx, Ly = $Ly, eltype = $T, usecuda = $usecuda.")
    nj_red = (Ly-1)*4 + 16 * Ly
    println("Layer 1/$Lx")
    k = 0
    reg = red_reg(T, Ly, Js[k+1:k+nj_red], hs[1:4Ly]; usecuda=usecuda)
    k += nj_red
    for j=1:Ly*4
        hs[j]!= 0 && (reg |> put(4*Ly, j=>Gh(T, hs[j])))
    end
    for i=2:Lx
        hk = (i-1)*Ly*8
        println("Layer $i/$Lx")
        # BLACK
        for j=1:Ly*4
            reg |> put(Ly*4, j=>G2(T, Js[k+1]))
            k += 1
        end
        for j=1:Ly*4 # `hs` interated in red->black order
            hs[j+hk+4Ly] != 0 && (reg |> put(4*Ly, j=>Gh(T, hs[j+hk+4Ly])))
        end

        # Contract with RED
        rr = red_reg(T, Ly, Js[k+1:k+nj_red], hs[hk+1:hk+4Ly]; usecuda=usecuda)
        k += nj_red
        reg.state .*= rr.state
    end
    if length(Js) !== k
        @warn "length of parameters is $(length(Js)), used $k"
    end
    sum(statevec(reg))
end

# NOTE: repeat with the definition in Viznet.
function sgbonds(lt::ChimeraLattice)
    bonds = Tuple{Int,Int}[]
    red_bond!(lt, 1, bonds)
    redblack_bond!(lt, 1, bonds)
    for i=2:lt.Nx
        # BLACK
        black_bond!(lt, i-1, bonds)
        # Contract with RED
        red_bond!(lt, i, bonds)
        redblack_bond!(lt, i, bonds)
    end
    return bonds
end

function red_bond!(lt, i::Int, bonds)
    CI = LinearIndices(lt |> size)
    for j=1:lt.Ny-1
        for k = 1:4
            push!(bonds, (CI[k, i, j], CI[k, i, j+1]))
        end
    end
    return bonds
end

function redblack_bond!(lt, i::Int, bonds)
    CI = LinearIndices(lt |> size)
    for j=1:lt.Ny
        for k2=5:8
            for k1=1:4
                push!(bonds, (CI[k2, i, j], CI[k1, i, j]))
            end
        end
    end
    return bonds
end

function black_bond!(lt, i::Int, bonds)
    CI = LinearIndices(lt |> size)
    for j=1:lt.Ny
        for k = 5:8
            push!(bonds, (CI[k,i,j], CI[k,i+1,j]))
        end
    end
    return bonds
end

# for beautiful printing.
function Base.getindex(lt::ChimeraLattice, k::Int, i::Int, j::Int)
    ki = (k-1) ÷ 4
    kj = mod1(k, 4)-1
    x = ((i-1)*(Viznet.gap_x(lt) + 2)+ki+0.5)*unit(lt)
    y = ((j-1)*(Viznet.gap_y(lt) + 4)+kj+0.5)*unit(lt)
    return (x+0.01*kj, y+0.02*ki)
end

# to index `h`.
function sgvertices(lt::ChimeraLattice)
    v = Int[]
    LI = LinearIndices(lt)
    for i=1:lt.Nx
        for ki = 1:2
            for j=1:lt.Ny
                for kj = 1:4
                    k = (ki-1)*4 + kj
                    push!(v, LI[k,i,j])
                end
            end
        end
    end
    return v
end

regsize(lt::ChimeraLattice) = lt.Ny*4
cachesize_A(lt::ChimeraLattice) = lt.Ny*4
cachesize_B(lt::ChimeraLattice) = lt.Nx-1
cachesize_largemem(lt::ChimeraLattice) = (lt.Nx-1)*lt.Ny*4

function _init_reg(::Type{T}, lt::ChimeraLattice, usecuda) where T
    _init_reg(T, lt.Ny*4, usecuda)
end
