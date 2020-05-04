using NiLang, NiLang.AD
using TropicalTensors
using BitBasis
using YaoArrayRegister: SDDiagonal, SDMatrix
using LinearAlgebra, Yao
using TupleTools
using StaticArrays

include("lib.jl")

@i function i_instruct!(state::StridedVector{<:Tropical},
        U0::AbstractMatrix{<:Tropical},
        locs::NTuple{M, Int},
        clocs::NTuple{C, Int},
        cvals::NTuple{C, Int}) where {C, M}
    @routine @invcheckoff begin
        nbit ← log2i(length(state))
        U ← (all(TupleTools.diff(locs).>0) ? U0 : reorder(U0, collect(locs)|>sortperm))
        MM ← size(U0, 1)
        locked_bits ← [clocs..., locs...]
        locked_vals ← [cvals..., zeros(Int, M)...]
        locs_raw ← [i+1 for i in itercontrol(nbit, setdiff(1:nbit, locs), zeros(Int, nbit-M))] |> YaoArrayRegister.staticize
        configs ← itercontrol(nbit, locked_bits, locked_vals)
    end
    loop_kernel(state, configs, U, locs_raw)
    ~@routine
end

@i function loop_kernel(state::StridedVector{<:Tropical}, configs, U::SDDiagonal, locs_raw)
    @invcheckoff @inbounds for i=1:length(configs)
        x ← configs[i]
        for i=1:length(U.diag)
            muleq(state[x+locs_raw[i]], U.diag[i])
        end
    end
end

@i function loop_kernel(state::StridedVector{T}, configs, U::AbstractMatrix, locs_raw) where T<:Tropical
    @invcheckoff begin
        branch_keeper ← zeros(Bool, size(U, 2))
        new_state ← ones(T, length(state))
        @inbounds for i=1:length(configs)
            x ← configs[i]
            #igemv!(view(new_state, x.+locs_raw), U, view(state, x.+locs_raw), branch_keeper)

			for l=1:size(U,1)
				el ← zero(T)
				@routine for k=1:size(U,2)
					muleq(U[l,k], state[x+locs_raw[k]])
					if (el.n < U[l,k].n, branch_keeper[k])
						FLIP(branch_keeper[k])
						NiLang.SWAP(el, U[l,k])
					end
				end
				muleq(new_state[x+locs_raw[l]], el)
				~@routine
			end
        end
        NiLang.SWAP(state, new_state)
        ipush!(new_state)
        new_state → zeros(T, length(state))
        branch_keeper → zeros(Bool, size(U, 2))
    end
end

@testset "new instr" begin
    g4 = Diagonal(spinglass_g4_tensor(1.5))
    reg = ArrayReg(randn(1<<12) .|> Tropical)
    s1 = i_instruct!(copy(vec(reg.state)), g4, (7, 3), (), ())[1]
    nreg = copy(reg) |> put(12, (7, 3)=>matblock(g4))
    @test statevec(nreg) ≈ s1

    g4 = randn(4, 4) .|> Tropical
    reg = ArrayReg(randn(1<<12) .|> Tropical)
    s1 = i_instruct!(copy(vec(reg.state)), g4, (7, 3), (), ())[1]
    nreg = copy(reg) |> put(12, (7, 3)=>matblock(g4))
    @test statevec(nreg) ≈ s1
end

using BenchmarkTools
g4 = Diagonal(MMatrix{4,4}(spinglass_g4_tensor(1.5)))
reg = ArrayReg(randn(1<<18) .|> Tropical)
@benchmark i_instruct!($(copy(vec(reg.state))), $g4, (3, 7), (), ())
@benchmark i_instruct!($(copy(GVar.(vec(reg.state)))), $(GVar(g4)), (3, 7), (), ())
@benchmark $(copy(reg)) |> $(put(18, (3, 7)=>matblock(g4)))

g4 = MMatrix{4,4}(Tropical.(randn(4,4)))
@benchmark (i_instruct!($(copy(vec(reg.state))), $g4, (3, 7), (), ()); empty!(NiLang.GLOBAL_STACK))
@benchmark (i_instruct!($(copy(GVar.(vec(reg.state)))), $(GVar(g4)), (3, 7), (), ()); empty!(NiLang.GLOBAL_STACK))
@benchmark copy($reg) |> $(put(18, (3, 7)=>matblock(g4)))
