using .CuYao
using .CuYao: CuArrays
using .CuArrays: CuArray

CuArrays.ones(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), one(Tropical{T}))
CuArrays.zeros(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), zero(Tropical{T}))
CuArrays.ones(::Type{CountingTropical{T}}, dims...) where T = fill!(CuArray{CountingTropical{T}}(undef, dims...), one(CountingTropical{T}))
CuArrays.zeros(::Type{CountingTropical{T}}, dims...) where T = fill!(CuArray{CountingTropical{T}}(undef, dims...), zero(CountingTropical{T}))

function _init_reg(::Type{T}, L::Int, usecuda::Val{:true}) where T
    ArrayReg(CuArrays.ones(T, 1<<L))
end

function _init_reg(::Type{T}, lt::MaskedSquareLattice, ::Val{:true}) where T
    nbit = size(lt, 2) + 2
    state = CuArrays.zeros(T, 1<<nbit)
    fill!(view(state,1:1<<(nbit-2)), one(T))
    ArrayReg(state)
end
