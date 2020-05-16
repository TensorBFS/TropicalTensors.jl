using ForwardDiff, Yao, LinearAlgebra

# define tropical numbers with the following property.
# x ⊕ y := max(x ,y)
# x ⊗ y := x + y
struct Tropical{T} <: Number
    n::T
end
Tropical{T}(x::Tropical{T}) where T = x

Base.:*(a::Tropical, b::Tropical) = Tropical(a.n + b.n)
Base.:+(a::Tropical, b::Tropical) = Tropical(max(a.n, b.n))
Base.zero(::Type{Tropical{T}}) where T<:Integer = Tropical(T(-999999))
Base.zero(::Type{Tropical{T}}) where T<:AbstractFloat = Tropical(typemin(T))

# define the "spinglass gates" that will be used in our "quantum" simulation.
# * `Gh` is the magnetic field term.
# * `G2` is the gate on parallel bond.
# * `G4` is the vertical two qubit diagonal "gate".
Gh(h::Real) = matblock(Diagonal(Tropical.([h, -h])))
G2(Jij::Real) = matblock(Tropical.([Jij -Jij; -Jij Jij]))
G4(Jij::Real) = matblock(Diagonal(Tropical.([Jij, -Jij, -Jij, Jij])))

function square_solve!(Lx::Int, Ly::Int, Js::AbstractVector{T}, hs::AbstractVector{T}; usecuda=false) where T
    Js, hs = copy(Js), copy(hs)
    reg = ArrayReg(Tropical.(zeros(T, 1<<Ly)))
    if usecuda
        reg = cu(reg)
    end
    for i=1:Lx
        println("Layer $i/$Lx")
        i!=1 && for j=1:Ly
            reg |> put(Ly, j=>G2(Js |> popfirst!))
        end
        for j=1:Ly
            reg |> put(Ly, j=>Gh(hs |> popfirst!))
        end
        for j=1:Ly-1
            reg |> put(Ly, (j,j+1)=>G4(Js |> popfirst!))
        end
    end
    sum(state(reg)).n
end

L = 15
Js = randn(2*L*(L-1))
hs = zeros(L^2)
gs = ForwardDiff.gradient(x->square_solve!(L, L,
    convert.(eltype(x), Js), x; usecuda=false), hs)

# use cuda
using CuYao
gs = ForwardDiff.gradient(x->square_solve!(L, L,
    convert.(eltype(x), Js), x; usecuda=true), hs)
