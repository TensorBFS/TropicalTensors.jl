using CuArrays

export Tropical

struct Tropical{T} <: Number
    n::T
end

Tropical{T}(x::Tropical{T}) where T = x

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
Base.zero(::Type{Tropical{T}}) where T = Tropical(T(-Inf))
Base.zero(::Type{Tropical{Int16}}) = Tropical(-32768)
Base.zero(::Type{Tropical{Int32}}) = Tropical(-2147483648)
Base.zero(::Type{Tropical{Int64}}) = Tropical(-9223372036854775808)
Base.zero(::Tropical{T}) where T = zero(Tropical{T})

Base.one(::Type{Tropical{T}}) where T = Tropical(zero(T))
Base.one(::Tropical{T}) where T = one(Tropical{T})

CuArrays.ones(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), one(Tropical{T}))
CuArrays.zeros(::Type{Tropical{T}}, dims...) where T = fill!(CuArray{Tropical{T}}(undef, dims...), zero(Tropical{T}))

#for OP in [:(Base.:>), :(Base.:<), :(Base.:==), :(Base.isapprox), :(Base.:>=), :(Base.:<=)]
for OP in [:>, :<, :(==), :>=, :<=]
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
