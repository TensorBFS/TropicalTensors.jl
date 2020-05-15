using Yao
using LuxurySparse
using LinearAlgebra
using TropicalYao
using Viznet

export solve, SquareLattice, ChimeraLattice

G2(::Type{T}, J) where T = matblock(spinglass_bond_tensor(T(J)) |> LuxurySparse.staticize)
G4(::Type{T}, J) where T = matblock(Diagonal(spinglass_g4_tensor(T(J))) |> LuxurySparse.staticize)
G16(::Type{T}, Js) where T = matblock(spinglass_g16_tensor(T.(Js)) |> LuxurySparse.staticize)

function _init_reg(::Type{T}, L::Int, usecuda::Val{:false}) where T
    reg = ArrayReg(ones(Tropical{T}, 1<<L))
end

include("square.jl")
include("chimera.jl")
