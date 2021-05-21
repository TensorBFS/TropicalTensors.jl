using .CuYao
using .CuYao: CUDA
using .CUDA: CuArray, @linearidx, GPUArrays
using LinearAlgebra

export togpu

CUDA.ones(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), one(Tropical{T}))
CUDA.zeros(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), zero(Tropical{T}))
CUDA.ones(::Type{CountingTropical{T,CT}}, dims...) where {T,CT} = fill!(CuArray{CountingTropical{T,CT}}(undef, dims...), one(CountingTropical{T,CT}))
CUDA.zeros(::Type{CountingTropical{T,CT}}, dims...) where {T,CT} = fill!(CuArray{CountingTropical{T,CT}}(undef, dims...), zero(CountingTropical{T,CT}))

function _init_reg(::Type{T}, L::Int, usecuda::Val{:true}) where T
    ArrayReg(CUDA.ones(T, 1<<L))
end

function _init_reg(::Type{T}, lt::MaskedSquareLattice, ::Val{:true}) where T
    nbit = size(lt, 2) + 2
    state = CUDA.ones(T, 1<<nbit)
    ArrayReg(state)
end

function togpu(tn::TensorNetwork)
    TensorNetwork(togpu.(tn.tensors))
end

function togpu(t::LabeledTensor)
    LabeledTensor(CuArray(t.array), t.labels)
end
