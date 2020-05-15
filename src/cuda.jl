using .CuYao: CuArrays
using CuArrays: CuArray

CuArrays.ones(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), one(Tropical{T}))
CuArrays.zeros(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), zero(Tropical{T}))

include("yao/cuda.jl")
