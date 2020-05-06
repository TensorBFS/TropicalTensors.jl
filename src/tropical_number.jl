export Tropical, TropicalF64, TropicalF32, TropicalF16

struct Tropical{T} <: Number
    n::T
    Tropical{T}(x) where T = new{T}(T(x))
    function Tropical(x::T) where T
        new{T}(x)
    end
end

const TropicalF64 = Tropical{Float64}
const TropicalF32 = Tropical{Float32}
const TropicalF16 = Tropical{Float16}

function Base.show(io::IO, inf::Tropical)
    print(io,"Tropical($(inf.n))")
end

function Base.show(io::IO, ::MIME"text/plain", inf::Tropical)
    Base.show(io, inf)
end

value(x::Tropical) = x.n

Base.:*(a::Tropical, b::Tropical) = Tropical(a.n + b.n)
function Base.:*(a::Tropical{<:Rational}, b::Tropical{<:Rational})
    if a == zero(a)
        a
    elseif b == zero(b)
        b
    else
        Tropical(a.n + b.n)
    end
end
Base.:+(a::Tropical, b::Tropical) = Tropical(max(a.n, b.n))
Base.zero(::Type{Tropical{T}}) where T<:Integer = Tropical(typemin(T)÷2)
Base.zero(::Type{Tropical{T}}) where T<:AbstractFloat = Tropical(typemin(T)/2)
Base.zero(::Tropical{T}) where T = zero(Tropical{T})

Base.one(::Type{Tropical{T}}) where T = Tropical(zero(T))
Base.one(::Tropical{T}) where T = one(Tropical{T})

#for OP in [:(Base.:>), :(Base.:<), :(Base.:==), :(Base.isapprox), :(Base.:>=), :(Base.:<=)]
for OP in [:>, :<, :(==), :>=, :<=, :isless]
    @eval Base.$OP(a::Tropical, b::Tropical) = $OP(a.n, b.n)
end

Base.isapprox(x::Tropical, y::Tropical; kwargs...) = isapprox(x.n, y.n; kwargs...)
Base.isapprox(x::AbstractArray{<:Tropical}, y::AbstractArray{<:Tropical}; kwargs...) = all(isapprox.(x, y; kwargs...))

# neighbor analysis
# degree == 2, from 1.2.2.2 and 1.2.1, we know the maximum branching is 2, <= 2^(1/3).
#

# TODO:
# 1. sparse tensor type to store possible configurations
# 2. merging two tensors, a flaw of this context un-aware algebra?
#    maybe it is also an advantage, it represents probability, gives a greedy version of algorithm.

# not working now, since gcd is not properly implemented for GPU, should be an easy fix.
#using CuArrays
#ein"ij,jk,ki->jk"(CuArray(T(1//2,1//2)), CuArray(T(1//2,1//2)), CuArray(T(1//2,1//2)))
