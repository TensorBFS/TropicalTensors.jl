export isolve, isolve_largemem

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
        for j = 1:length(reg.state)
            @inbounds B_STACK[j] *= identity(reg.state[j])
        end

        # clean up `NiLang.GLOBAL_STACK`
        ~@routine
        # but wait, `k` should not be uncomputed
        k += identity(Ly)

        # store the stored state
        for j=1:length(reg.state)
            @inbounds NiLang.SWAP(reg.state[j], B_STACK[j])
        end

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

function opt_config(sg::Spinglass{LT,T}) where {LT,T}
    reg = ArrayReg(ones(Tropical{T}, 1<<sg.lattice.Ny))
    A = stack4reg(reg, sg.lattice.Ny)
    B = stack4reg(reg, sg.lattice.Nx-1)
    eng, sg, reg, A, B = isolve(T(0.0), sg, reg, A, B)
    sgg = Spinglass(sg.lattice, GVar.(sg.Js, zero(sg.Js)), GVar.(sg.hs, zero(sg.hs)))
    gres = (~isolve)(GVar(eng, T(1)), sgg, GVar(reg), GVar(A), GVar(B))
    empty!(NiLang.GLOBAL_STACK)
    return eng, grad.(gres[2].Js)
end
