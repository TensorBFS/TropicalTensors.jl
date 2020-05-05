using YaoBlocks
using StaticArrays

export tropicalblock, TropicalMatrixBlock
export apply_G4!, apply_G2!, apply_G16!

"""
    TropicalMatrixBlock{N, MT} <: PrimitiveBlock{N}

Tropical matrix gate wraps a tropical matrix operator to Yao blocks.
`N` is the hilbert dimension of current block.
"""
struct TropicalMatrixBlock{N,T,MT<:AbstractMatrix{T}} <: PrimitiveBlock{N}
    mat::MT
end

YaoBlocks.PutBlock{N,C,GT}(g::GT, locs) where {N, C, GT} = PutBlock{N}(g, locs)

function TropicalMatrixBlock{N}(m::MT) where {N,T,MT<:AbstractMatrix{T}}
    (1 << N, 1 << N) == size(m) ||
    throw(DimensionMismatch("expect a $(1<<N) x $(1<<N) matrix, got $(size(m))"))
    return TropicalMatrixBlock{N,T,MT}(m)
end

TropicalMatrixBlock(m::AbstractMatrix) = TropicalMatrixBlock{log2i(size(m,1))}(m)
(f::Type{Inv{TropicalMatrixBlock}})(m::TropicalMatrixBlock) = m.mat

tropicalblock(m::AbstractMatrix{<:Tropical}) = TropicalMatrixBlock(m)
mat(A::TropicalMatrixBlock) = A.mat

function mat(::Type{T}, A::TropicalMatrixBlock) where {T}
    if eltype(A.mat) == T
        return A.mat
    else
        # this errors before, but since we allow one to specify T in mat
        # this should be allowed but with a suggestion
        @warn "converting $(eltype(A.mat)) to eltype $T, consider create another matblock with eltype $T"
        return copyto!(similar(A.mat, T), A.mat)
    end
end

Base.:(==)(A::TropicalMatrixBlock, B::TropicalMatrixBlock) = A.mat == B.mat
Base.copy(A::TropicalMatrixBlock) = TropicalMatrixBlock(copy(A.mat))
Base.adjoint(x::TropicalMatrixBlock) = Daggered(x)

function NiLang.chfield(a::AbstractArray, ::typeof(vec), val)
    reshape(val, size(a)...)
end

@i function YaoBlocks.apply!(reg::ArrayReg{1}, pb::PutBlock{N,C,<:TropicalMatrixBlock}) where {N, C}
    i_instruct!(vec(reg.state), pb.content.mat, pb.locs, (), ())
end
YaoBlocks.apply!(reg::ArrayReg{B}, b::TropicalMatrixBlock) where B = throw(NotImplementedError(:apply!, typeof((reg, b))))

@i function apply_G2!(reg::ArrayReg{1,T}, i::Int, J::Real) where T<:Tropical
    @routine begin
        nbit ← nqubits(reg)
        blk ← put(nbit, i=>tropicalblock(MMatrix{2,2}(ones(T, 2, 2))))
        spinglass_bond_tensor!(blk.content.mat, J)
    end
    apply!(reg, blk)
    ~@routine
end

@i function apply_G4!(reg::ArrayReg{1,T}, i::NTuple{2,Int}, J::Real) where T<:Tropical
    @routine begin
        nbit ← nqubits(reg)
        blk ← put(nbit, i=>tropicalblock(Diagonal(MVector{4}(ones(T, 4)))))
        spinglass_g4_tensor!(blk.content.mat, J)
    end
    apply!(reg, blk)
    ~@routine
end

@i function apply_G16!(reg::ArrayReg{1,T}, i::NTuple{4,Int}, Js::AbstractVector{<:Real}) where T<:Tropical
    @routine begin
        nbit ← nqubits(reg)
        blk ← put(nbit, i=>tropicalblock(ones(T, 16, 16)))
        spinglass_g16_tensor!(blk.content.mat, Js)
    end
    apply!(reg, blk)
    ~@routine
end
