using .CuYao
using .CuYao: CuArrays
using .CuArrays: CuArray

CuArrays.ones(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), one(Tropical{T}))
CuArrays.zeros(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), zero(Tropical{T}))

function _init_reg(::Type{T}, L::Int, usecuda::Val{:true}) where T
    ArrayReg(CuArrays.ones(Tropical{T}, 1<<L))
end

function _init_reg(::Type{T}, lt::MaskedSquareLattice, ::Val{:true}) where T
    nbit = size(lt, 2) + 2
    state = CuArrays.zeros(Tropical{T}, 1<<nbit)
    fill!(view(state,1:1<<(nbit-2)), one(Tropical{T}))
    ArrayReg(state)
end
