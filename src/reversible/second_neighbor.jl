function _c(lt, a, b)
    LI = LinearIndices(lt)
    isconnected(lt, LI[a...], LI[b...])
end

@i function isolve_largemem(out!::T, sg::AbstractSpinglass{MaskedSquareLattice}, reg::ArrayReg{1,TT}, REG_STACK; usecuda=false) where {T,TT<:Tropical{T}}
    @invcheckoff begin
    lt ← sg.lattice
    Lx ← size(lt, 1)
    Ly ← size(lt, 2)
    nbit ← Ly + 2
    Js ← sg.Js
    hs ← sg.hs
    k ← 0
    l ← 0
    for i=1:Lx
        @safe println("Layer $i/$Lx")
        for j=1:Ly-1
            if (_c(lt, (i,j), (i,j+1)), ~)
                INC(k)
                apply_G4!(reg, (j, j+1), Js[k], REG_STACK)
            end
        end
        for j=1:Ly
            if (lt.mask[i,j], ~)
                INC(l)
                apply_Gh!(reg, j, hs[l], REG_STACK)
            end
        end
        if (i!=Lx, ~)
            for j=1:Ly
                # store the information in qubit `j` to ancilla `nbit-j%2`
                if (j!=Ly && _c(lt, (i,j), (i+1,j+1)), ~)
                    apply_Gcp!(reg, (j,nbit-j%2), REG_STACK)
                end
                # interact with j-1 th qubit (a)
                if (j!=1 && _c(lt, (i+1,j-1), (i,j)), ~)
                    INC(k)
                    apply_G4!(reg, (j-1, j), Js[k], REG_STACK)
                end
                # onsite term (b)
                if (_c(lt, (i,j), (i+1,j)), ~)
                    INC(k)
                    apply_G2!(reg, j, Js[k], REG_STACK)
                else
                    apply_Gcut!(reg, j, REG_STACK)
                end
                if (j!=1 && _c(lt, (i,j-1), (i+1,j)), ~)
                    INC(k)
                    # interact with cached j-1 th qubit (c)
                    apply_G4!(reg, (j,nbit-(j-1)%2), Js[k], REG_STACK)
                    # erease the information in previous ancilla `nbit-(j-1)%2`
                    apply_Gcut!(reg, nbit-(j-1)%2, REG_STACK)
                end
            end
        end
    end
    summed ← one(TT)
    isum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    l → length(hs)
    k → length(Js)
    summed → one(TT)
    end
end

function cachesize_largemem(lt::MaskedSquareLattice)
    Lx = size(lt, 1)
    Ly = size(lt, 2)
    ncache = 0
    for i=1:Lx
        i!=Lx && for j=1:Ly
            if j!=Ly && _c(lt, (i,j), (i+1,j+1))
                ncache += 1
            end
            ncache += 1
            if j!=1 && _c(lt, (i,j-1), (i+1,j))
                ncache += 1
            end
        end
    end
    return ncache
end

@i function isolve(out!::T, sg::Spinglass{MaskedSquareLattice,T}, reg::ArrayReg{1,TT}, A_STACK, B_STACK; usecuda=false) where {T,TT<:Tropical{T}}
    lt ← sg.lattice
    Lx ← size(lt, 1)
    Ly ← size(lt, 2)
    nbit ← Ly + 2
    Js ← sg.Js
    hs ← sg.hs
    k ← 0
    l ← 0
    for i=1:Lx
        @safe println("Layer $i/$Lx")
        for j=1:Ly-1
            if (_c(lt, (i,j), (i,j+1)), ~)
                INC(k)
                apply_G4!(reg, (j, j+1), Js[k], A_STACK)
            end
        end
        for j=1:Ly
            if (lt.mask[i,j], ~)
                INC(l)
                apply_Gh!(reg, j, hs[l], A_STACK)
            end
        end
        if (i!=Lx, ~)
            @routine begin
                for j=1:Ly
                    # store the information in qubit `j` to ancilla `nbit-j%2`
                    if (j!=Ly && _c(lt, (i,j), (i+1,j+1)), ~)
                        apply_Gcp!(reg, (j,nbit-j%2), A_STACK)
                    end
                    # interact with j-1 th qubit (a)
                    if (j!=1 && _c(lt, (i+1,j-1), (i,j)), ~)
                        INC(k)
                        apply_G4!(reg, (j-1, j), Js[k], A_STACK)
                    end
                    # onsite term (b)
                    if (_c(lt, (i,j), (i+1,j)), ~)
                        INC(k)
                        apply_G2!(reg, j, Js[k], A_STACK)
                    else
                        apply_Gcut!(reg, j, A_STACK)
                    end
                    if (j!=1 && _c(lt, (i,j-1), (i+1,j)), ~)
                        INC(k)
                        # interact with cached j-1 th qubit (c)
                        apply_G4!(reg, (j,nbit-(j-1)%2), Js[k], A_STACK)
                        # erease the information in previous ancilla `nbit-(j-1)%2`
                        apply_Gcut!(reg, nbit-(j-1)%2, A_STACK)
                    end
                end
            end
            incstack!(B_STACK)
            store_state!(B_STACK, reg.state)
            ~@routine
            # but wait, `k` should not be uncomputed
            for j=1:Ly
                if (j!=1 && _c(lt, (i+1,j-1), (i,j)), ~)
                    INC(k)
                end
                # onsite term (b)
                if (_c(lt, (i,j), (i+1,j)), ~)
                    INC(k)
                end
                if (j!=1 && _c(lt, (i,j-1), (i+1,j)), ~)
                    INC(k)
                end
            end
            # restore state
            swap_state!(B_STACK, reg.state)
        end
    end
    summed ← one(TT)
    isum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    l → length(hs)
    k → length(Js)
    summed → one(TT)
end

cachesize_A(lt::MaskedSquareLattice) = size(lt,2)*3-1
cachesize_B(lt::MaskedSquareLattice) = size(lt,1)-1
