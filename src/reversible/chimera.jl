@i function red_reg(reg::ArrayReg{B,T}, Ly::Int, Js, hs, REG_STACK) where {B,T}
    k ← 0
    for j=1:Ly*4
        apply_Gh!(reg, j, hs[j], REG_STACK)
    end
    @invcheckoff for i=1:Ly-1
        for j = 1:4
            apply_Gvb!(reg, (@const (4*(i-1)+j,4*i+j)), Js[k+1], REG_STACK)
            k += 1
        end
    end
    @invcheckoff for i=1:Ly
        apply_G16!(reg, (@const (4i-3,4i-2,4i-1,4i)), Js[k+1:k+16], REG_STACK)
        k += 16
    end
    k → length(Js)
end

@i function isolve_largemem(out!::T, sg::Spinglass{ChimeraLattice, T}, reg::ArrayReg{B,TT}, REG_STACK) where {B,T,TT<:Tropical{T}}
    @invcheckoff begin
    @routine begin
        Lx ← sg.lattice.Nx
        Ly ← sg.lattice.Ny
        Js ← sg.Js
        hs ← sg.hs
        nj_red ← (Ly-1)*4 + 16 * Ly
    end
    k ← 0
    @safe println("Running Chimera, Lx = $Lx, Ly = $Ly, eltype = $T.")
    @safe println("Layer 1/$Lx")
    red_reg(reg, Ly, Js[k+1:k+nj_red], hs[1:4Ly], REG_STACK)
    k += nj_red
    for j=1:Ly*4 # `hs` interated in red->black order
        apply_Gh!(reg, j, hs[j+4Ly], REG_STACK)
    end

    for i=2:Lx
        @routine hk ← (i-1)*Ly*8
        @safe println("Layer $i/$Lx")
        # BLACK
        for j=1:Ly*4
            apply_Ghb!(reg, j, Js[k+1], REG_STACK)
            k += 1
        end
        for j=1:Ly*4 # `hs` interated in red->black order
            apply_Gh!(reg, j, hs[j+hk+4Ly], REG_STACK)
        end

        # Contract with RED
        @routine rr ← _init_reg(Tropical{T}, Ly*4, Val(false))
        red_reg(rr, Ly, Js[k+1:k+nj_red], hs[hk+1:hk+4Ly], REG_STACK)
        reg.state .*= rr.state
        incstack!(REG_STACK)
        swap_state!(REG_STACK, rr.state)
        k += nj_red
        ~@routine
        ~@routine
    end
    summed ← one(TT)
    i_sum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(Js)
    summed → one(TT)
    ~@routine
    end
end

@i function isolve(out!::T, sg::AbstractSpinglass{ChimeraLattice, T}, reg::ArrayReg{B,TT}, A_STACK, B_STACK) where {B,T,TT<:Tropical{T}}
    @invcheckoff begin
    @routine begin
        Lx ← sg.lattice.Nx
        Ly ← sg.lattice.Ny
        Js ← sg.Js
        hs ← sg.hs
        nj_red ← (Ly-1)*4 + 16 * Ly
    end
    k ← 0
    @safe println("Running Chimera, Lx = $Lx, Ly = $Ly, eltype = $T.")
    @safe println("Layer 1/$Lx")

    # initial reg
    @routine begin
        rr ← _init_reg(Tropical{T}, Ly*4, Val(false))
        red_reg(rr, Ly, Js[k+1:k+nj_red], hs[1:4Ly], A_STACK)
    end
    reg.state .*= rr.state
    ~@routine
    k += nj_red

    for j=1:Ly*4 # `hs` interated in red->black order
        apply_Gh!(reg, j, hs[j+4Ly], A_STACK)
    end

    for i=2:Lx
        @safe println("Layer $i/$Lx")
        hk ← (i-1)*Ly*8
        # BLACK
        @routine begin
            for j=1:Ly*4
                apply_Ghb!(reg, j, Js[k+1], A_STACK)
                k += 1
            end
        end
        incstack!(B_STACK)
        store_state!(B_STACK, reg.state)
        ~@routine
        k += 4*Ly
        swap_state!(B_STACK, reg.state)

        for j=1:Ly*4 # `hs` interated in red->black order
            apply_Gh!(reg, j, hs[j+hk+4Ly], A_STACK)
        end

        # RED
        @routine begin
            rr ← _init_reg(Tropical{T}, Ly*4, Val(false))
            red_reg(rr, Ly, Js[k+1:k+nj_red], hs[hk+1:hk+4Ly], A_STACK)
        end
        reg.state .*= rr.state
        ~@routine
        k += nj_red
        hk → (i-1)*Ly*8
    end
    summed ← one(TT)
    i_sum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(Js)
    summed → one(TT)
    ~@routine
    end
end

cachesize_A(lt::ChimeraLattice) = lt.Ny*4
cachesize_B(lt::ChimeraLattice) = lt.Nx-1
cachesize_largemem(lt::ChimeraLattice) = (lt.Nx-1)*lt.Ny*4
