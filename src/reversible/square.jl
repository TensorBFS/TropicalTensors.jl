@i function isolve(out!::T, sg::AbstractSpinglass{SquareLattice, T}, reg::ArrayReg{B,TT}, A_STACK, B_STACK) where {B,T,TT<:Tropical{T}}
    @invcheckoff begin
    Lx ← sg.lattice.Nx
    Ly ← sg.lattice.Ny
    Js ← sg.Js
    hs ← sg.hs
    @safe println("Layer 1/$Lx, stack size: $(A_STACK.top) & $(B_STACK.top)")
    k ← 0
    for j=1:Ly
        apply_Gh!(reg, j, hs[j], A_STACK)
    end
    for j=1:Ly-1
        k += 1
        apply_Gvb!(reg, (@const (j, j+1)), Js[k], A_STACK)
    end
    for i=2:Lx
        @safe println("Layer $i/$Lx, stack size: $(A_STACK.top) & $(B_STACK.top)")
        @routine begin
            for j=1:Ly
                k += 1
                apply_Ghb!(reg, j, Js[k], A_STACK)
            end
        end
        # bookup current state and loss
        incstack!(B_STACK)
        store_state!(B_STACK, reg.state)

        # clean up `NiLang.GLOBAL_STACK`
        ~@routine
        # but wait, `k` should not be uncomputed
        k += Ly

        # restore the state
        swap_state!(B_STACK, reg.state)

        for j=1:Ly
            apply_Gh!(reg, j, hs[(i-1)*Ly+j], A_STACK)
        end
        for j=1:Ly-1
            k += 1
            apply_Gvb!(reg, (@const (j, j+1)), Js[k], A_STACK)
        end
    end
    summed ← one(TT)
    i_sum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(Js)
    summed → one(TT)
    end
end

@i function swap_state!(B_STACK, state)
    @invcheckoff for j=1:length(state)
        @inbounds NiLang.SWAP(state[j], B_STACK[j])
    end
end

@i function store_state!(B_STACK, state)
    @invcheckoff for j = 1:length(state)
        @inbounds TropicalYao.Reversible.unsafe_store!(B_STACK[j], state[j])
    end
end

@i function isolve_largemem(out!, sg::Spinglass{SquareLattice, T}, reg::ArrayReg{B,TT}, REG_STACK) where {B,T,TT<:Tropical{T}}
    @invcheckoff begin
    Lx ← sg.lattice.Nx
    Ly ← sg.lattice.Ny
    Js ← sg.Js
    hs ← sg.hs
    @safe println("Layer 1/$Lx")
    k ← 0
    for j=1:Ly
        apply_Gh!(reg, j, hs[j], REG_STACK)
    end
    for j=1:Ly-1
        k += 1
        apply_Gvb!(reg, (@const (j, j+1)), Js[k], REG_STACK)
    end
    for i=2:Lx
        @safe println("Layer $i/$Lx")
        for j=1:Ly
            k += 1
            apply_Ghb!(reg, j, Js[k], REG_STACK)
        end
        for j=1:Ly
            apply_Gh!(reg, j, hs[(i-1)*Ly+j], REG_STACK)
        end
        for j=1:Ly-1
            k += 1
            apply_Gvb!(reg, (@const (j, j+1)), Js[k], REG_STACK)
        end
    end
    summed ← one(TT)
    i_sum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(Js)
    summed → one(TT)
    end
end

cachesize_A(lt::SquareLattice) = lt.Ny
cachesize_B(lt::SquareLattice) = lt.Nx-1
cachesize_largemem(lt::SquareLattice) = (lt.Nx-1) * lt.Ny
