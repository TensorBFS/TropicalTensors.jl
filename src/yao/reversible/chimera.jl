@i function red_reg(reg::ArrayReg{B,T}, Ly::Int, Js, REG_STACK) where {B,T}
    k ← 0
    @invcheckoff for i=1:Ly-1
        for j = 1:4
            apply_G4!(reg, (4*(i-1)+j,4*i+j), Js[k+1], REG_STACK)
            k += identity(1)
        end
    end
    @invcheckoff for i=1:Ly
        apply_G16!(reg, (4i-3,4i-2,4i-1,4i), Js[k+1:k+16], REG_STACK)
        k += identity(16)
    end
    k → length(Js)
end

@i function isolve_largemem(out!::T, sg::Spinglass{<:ChimeraLattice, T}, reg::ArrayReg{B,TT}, REG_STACK) where {B,T,TT<:Tropical{T}}
    @invcheckoff begin
    Lx ← sg.lattice.Nx
    Ly ← sg.lattice.Ny
    J ← sg.Js
    nj_red ← (Ly-1)*4 + 16 * Ly
    k ← 0
    @safe println("Running Chimera, Lx = $Lx, Ly = $Ly, eltype = $T.")
    @safe println("Layer 1/$Lx")
    red_reg(reg, Ly, J[k+1:k+nj_red], REG_STACK)
    k += identity(nj_red)

    for j=2:Lx
        @safe println("Layer $j/$Lx")
        # BLACK
        for i=1:Ly*4
            apply_G2!(reg, i, J[k+1], REG_STACK)
            k += identity(1)
        end

        # Contract with RED
        rr ← _init_reg(T, Ly*4, Val(false))
        red_reg(rr, Ly, J[k+1:k+nj_red], REG_STACK)
        reg.state .*= identity.(rr.state)
        incstack!(REG_STACK)
        swap_state!(REG_STACK, rr.state)
        k += identity(nj_red)
    end
    summed ← one(TT)
    isum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(J)
    summed → one(TT)
    end
end

@i function isolve(out!::T, sg::Spinglass{<:ChimeraLattice, T}, reg::ArrayReg{B,TT}, A_STACK, B_STACK) where {B,T,TT<:Tropical{T}}
    @invcheckoff begin
    Lx ← sg.lattice.Nx
    Ly ← sg.lattice.Ny
    J ← sg.Js
    nj_red ← (Ly-1)*4 + 16 * Ly
    k ← 0
    @safe println("Running Chimera, Lx = $Lx, Ly = $Ly, eltype = $T.")
    @safe println("Layer 1/$Lx")

    # initial reg
    @routine begin
        rr ← _init_reg(T, Ly*4, Val(false))
        red_reg(rr, Ly, J[k+1:k+nj_red], A_STACK)
    end
    reg.state .*= identity.(rr.state)
    ~@routine

    k += identity(nj_red)

    for j=2:Lx
        @safe println("Layer $j/$Lx")
        # BLACK
        @routine begin
            for i=1:Ly*4
                apply_G2!(reg, i, J[k+1], A_STACK)
                k += identity(1)
            end
        end
        incstack!(B_STACK)
        store_state!(B_STACK, reg.state)
        ~@routine
        k += 4*Ly
        swap_state!(B_STACK, reg.state)

        # RED
        @routine begin
            rr ← _init_reg(T, Ly*4, Val(false))
            red_reg(rr, Ly, J[k+1:k+nj_red], A_STACK)
        end
        reg.state .*= identity.(rr.state)
        ~@routine
        k += identity(nj_red)
    end
    summed ← one(TT)
    isum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(J)
    summed → one(TT)
    end
end
