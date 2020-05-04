using NiLang
using TropicalTensors
using BitBasis

include("lib.jl")

function YaoBase.instruct!(
    state::AbstractVecOrMat{T1},
    operator::AbstractMatrix{T2},
    locs::NTuple{M,Int},
    control_locs::NTuple{C,Int} = (),
    control_bits::NTuple{C,Int} = (),
) where {T1,T2,M,C}

    if T2 != T1
        @warn "Element Type Mismatch: register $(T1), operator $(T2). Converting operator to match, this may cause performance issue"
        operator = copyto!(similar(operator, T1), operator)
    end
    operator = sort_unitary(operator, locs)
    locs_raw, ic = _prepare_instruct(state, operator, locs, control_locs, control_bits)
    return _instruct!(state, autostatic(operator), locs_raw, ic)
end

function _instruct!(
    state::AbstractVecOrMat{T},
    U::AbstractMatrix{T},
    locs_raw::SVector,
    ic::IterControl,
) where {T}
    controldo(ic) do i
        @inbounds unrows!(state, locs_raw .+ i, U)
    end
    return state
end

@i function instruct!(state::AbstractVector{<:Tropical},
        U0::AbstractMatrix{<:Tropical},
        locs::NTuple{M, Int},
        clocs::NTuple{C, Int},
        cvals::NTuple{C, Int}) where {C, M}
    nbit ← log2i(length(state))
    U ← (all(TupleTools.diff(locs).>0) ? U0 : reorder(U0, collect(locs)|>sortperm))
    MM ← size(U0, 1)
    locked_bits ← [clocs..., locs...]
    locked_vals ← [cvals..., zeros(Int, M)...]
    locs_raw ← [i+1 for i in itercontrol(nbit, setdiff(1:nbit, locs), zeros(Int, nbit-M))] |> YaoArrayRegister.staticize
    configs ← itercontrol(nbit, locked_bits, locked_vals)

    for i=1:length(configs)
        x ← @inbounds configs[i]
        unrows!(piecewise(reg.state, inds), x .+ locs_raw, U)
        igemv!(out[inds], U, state[inds], branch_keeper)
    end
end

@i function i_instruct!(state::AbstractVector{<:Tropical},
        U0::SDDiagonal{<:Tropical},
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

@i function loop_kernel(state, configs, U, locs_raw)
    @invcheckoff @inbounds for i=1:length(configs)
        x ← configs[i]
        for i=1:length(U.diag)
            muleq(state[x+locs_raw[i]], U.diag[i])
        end
    end
end

using YaoArrayRegister: SDDiagonal
using LinearAlgebra, Yao
using TupleTools

using NiLang, NiLang.AD
NiLang.AD.GVar(x::Tropical) = Tropical(GVar(x.n))
(_::Type{Inv{GVar}})(x::Tropical{<:GVar}) = Tropical((~GVar)(x.n))
NiLang.AD.GVar(x::ArrayReg{B}) where B = ArrayReg{B}(GVar(ArrayReg.state))

@testset "GVar" begin
    x = Tropical(0.4)
    @test GVar(x) isa Tropical{<:GVar}
    @test (~GVar)(GVar(Tropical(0.4))) == x
end

@testset "new instr" begin
    g4 = Diagonal(spinglass_g4_tensor(1.5))
    reg = ArrayReg(randn(1<<12) .|> Tropical)
    s1 = i_instruct!(copy(vec(reg.state)), g4, (3, 7), (), ())[1]
    nreg = copy(reg) |> put(12, (3, 7)=>matblock(g4))
    @test statevec(nreg) ≈ s1
end

g4 = Diagonal(MMatrix{4,4}(spinglass_g4_tensor(1.5)))
reg = ArrayReg(randn(1<<18) .|> Tropical)
@benchmark i_instruct!($(copy(vec(reg.state))), $g4, (3, 7), (), ())
@benchmark i_instruct!($(copy(GVar.(vec(reg.state)))), $(GVar(g4)), (3, 7), (), ())
@benchmark $(copy(reg)) |> $(put(18, (3, 7)=>matblock(g4)))
