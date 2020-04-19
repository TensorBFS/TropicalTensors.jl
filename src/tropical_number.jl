export Tropical

const TROPICAL_ZERO = -999999

struct Tropical{T} <: Number
    n::T
end

function Base.show(io::IO, inf::Tropical)
    print(io,"Tropical($(inf.n))")
end

function Base.show(io::IO, ::MIME"text/plain", inf::Tropical)
    Base.show(io, inf)
end

Base.:*(a::Tropical, b::Tropical) = Tropical(a.n + b.n)
Base.:+(a::Tropical, b::Tropical) = Tropical(max(a.n, b.n))
Base.zero(::Type{Tropical{T}}) where T = Tropical(T(TROPICAL_ZERO))
Base.zero(::Tropical{T}) where T = zero(Tropical{T})
#for OP in [:(Base.:>), :(Base.:<), :(Base.:==), :(Base.isapprox), :(Base.:>=), :(Base.:<=)]
for OP in [:>, :<, :(==), :(isapprox), :>=, :<=]
    @eval Base.$OP(a::Tropical, b::Tropical) = $OP(a.n, b.n)
end

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
