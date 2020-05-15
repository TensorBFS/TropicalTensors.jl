function red_reg(::Type{T}, Ly::Int, Js; usecuda=false) where T
    reg = _init_reg(T, Ly*4, Val(usecuda))
    k = 0
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
    J = sg.Js
    println("Running Chimera, Lx = $Lx, Ly = $Ly, eltype = $T, usecuda = $usecuda.")
    nj_red = (Ly-1)*4 + 16 * Ly
    println("Layer 1/$Lx")
    k = 0
    reg = red_reg(T, Ly, J[k+1:k+nj_red]; usecuda=usecuda)
    k += nj_red

    for j=2:Lx
        println("Layer $j/$Lx")
        # BLACK
        for i=1:Ly*4
            reg |> put(Ly*4, i=>G2(T, J[k+1]))
            k += 1
        end

        # Contract with RED
        rr = red_reg(T, Ly, J[k+1:k+nj_red]; usecuda=usecuda)
        k += nj_red
        reg.state .*= rr.state
    end
    if length(J) !== k
        @warn "length of parameters is $(length(J)), used $k"
    end
    sum(statevec(reg))
end

sgbonds(lt::ChimeraLattice) = bonds(lt)
