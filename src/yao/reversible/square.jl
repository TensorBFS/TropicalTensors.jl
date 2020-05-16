@i function isolve(out!, sg::Spinglass{<:SquareLattice, T}, reg::ArrayReg{B,TT}, A_STACK, B_STACK) where {B,T,TT<:Tropical{T}}
    @invcheckoff begin
    Lx ← sg.lattice.Nx
    Ly ← sg.lattice.Ny
    J ← sg.Js
    @safe println("Layer 1/$Lx, stack size: $(A_STACK.top) & $(B_STACK.top)")
    k ← 0
    for j=1:Ly-1
        k += identity(1)
        apply_G4!(reg, (j, j+1), J[k], A_STACK)
    end
    for i=2:Lx
        @safe println("Layer $i/$Lx, stack size: $(A_STACK.top) & $(B_STACK.top)")
        @routine begin
            for j=1:Ly
                k += identity(1)
                apply_G2!(reg, j, J[k], A_STACK)
            end
        end
        # bookup current state and loss
        incstack!(B_STACK)
        store_state!(B_STACK, reg.state)

        # clean up `NiLang.GLOBAL_STACK`
        ~@routine
        # but wait, `k` should not be uncomputed
        k += identity(Ly)

        # store the stored state
        swap_state!(B_STACK, reg.state)

        for j=1:Ly-1
            k += identity(1)
            apply_G4!(reg, (j, j+1), J[k], A_STACK)
        end
    end
    summed ← one(TT)
    isum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(J)
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
        @inbounds B_STACK[j] *= identity(state[j])
    end
end

@i function isolve_largemem(out!, sg::Spinglass{<:SquareLattice, T}, reg::ArrayReg{B,TT}, REG_STACK) where {B,T,TT<:Tropical{T}}
    @invcheckoff begin
    Lx ← sg.lattice.Nx
    Ly ← sg.lattice.Ny
    J ← sg.Js
    @safe println("Layer 1/$Lx")
    k ← 0
    for j=1:Ly-1
        k += identity(1)
        apply_G4!(reg, (j, j+1), J[k], REG_STACK)
    end
    for i=2:Lx
        @safe println("Layer $i/$Lx")
        for j=1:Ly
            k += identity(1)
            apply_G2!(reg, j, J[k], REG_STACK)
        end
        for j=1:Ly-1
            k += identity(1)
            apply_G4!(reg, (j, j+1), J[k], REG_STACK)
        end
    end
    summed ← one(TT)
    isum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(J)
    summed → one(TT)
    end
end
