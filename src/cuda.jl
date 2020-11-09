using .CuYao
using .CuYao: CUDA
using .CUDA: CuArray, @linearidx, GPUArrays

export togpu

CUDA.ones(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), one(Tropical{T}))
CUDA.zeros(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), zero(Tropical{T}))
CUDA.ones(::Type{CountingTropical{T}}, dims...) where T = fill!(CuArray{CountingTropical{T}}(undef, dims...), one(CountingTropical{T}))
CUDA.zeros(::Type{CountingTropical{T}}, dims...) where T = fill!(CuArray{CountingTropical{T}}(undef, dims...), zero(CountingTropical{T}))
Base.:(*)(a::CountingTropical, b::Bool) = b ? a : zero(a)
Base.:(*)(b::Bool, a::CountingTropical) = b ? a : zero(a)

function _init_reg(::Type{T}, L::Int, usecuda::Val{:true}) where T
    ArrayReg(CUDA.ones(T, 1<<L))
end

function _init_reg(::Type{T}, lt::MaskedSquareLattice, ::Val{:true}) where T
    nbit = size(lt, 2) + 2
    state = CUDA.zeros(T, 1<<nbit)
    fill!(view(state,1:1<<(nbit-2)), one(T))
    ArrayReg(state)
end

function togpu(tn::TensorNetwork)
    TensorNetwork(togpu.(tn.tensors); metas=tn.metas)
end

function togpu(t::LabeledTensor)
    LabeledTensor(CuArray(t.array), t.labels)
end

function GPUArrays.genperm(I::NTuple{N}, perm::NTuple{N}) where N
    ntuple(d-> (@inbounds return I[perm[d]]), Val(N))
end

function LinearAlgebra.permutedims!(dest::GPUArrays.AbstractGPUArray, src::GPUArrays.AbstractGPUArray, perm)
    perm isa Tuple || (perm = Tuple(perm))
    size_dest = size(dest)
    size_src = size(src)
    CUDA.gpu_call(vec(dest), vec(src), perm; name="permutedims!") do ctx, dest, src, perm
        i = @linearidx src
        I = l2c(size_src, i)
        @inbounds dest[c2l(size_dest, GPUArrays.genperm(I, perm))] = src[i]
        return
    end
    return reshape(dest, size(dest))
end

using Base.Cartesian
@generated function c2l(size::NTuple{N, Int}, c::NTuple{N,Int}) where N
    quote
        res = c[1]
        stride = size[1]
        @nexprs $(N-1) i->begin
            res += (c[i+1]-1) * stride
            stride *= size[i+1]
        end
        return res
    end
end

@generated function l2c(size::NTuple{N, Int}, l::Int) where N
    quote
        l -= 1
        @nexprs $(N-1) i->begin
            s_i = l % size[i] + 1
            l = l ÷ size[i]
        end
        $(Expr(:tuple, [Symbol(:s_, i) for i=1:N-1]..., :(l+1)))
    end
end