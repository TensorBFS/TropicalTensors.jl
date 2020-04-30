using Test
using Yao
using CuYao
using LinearAlgebra
using CuArrays
CuArrays.allowscalar(false)
using TropicalTensors
using LuxurySparse

_get_J(::Val{:ferro}) = 1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()

function _spinglass_yao(reg::ArrayReg{B,Tropical{T}}, L::Int, jtype::Val) where {B,T}
    G2 = matblock(spinglass_bond_tensor(T(1.0)) |> LuxurySparse.staticize)
    G4 = matblock(Diagonal(spinglass_g4_tensor(T(1.0))) |> LuxurySparse.staticize)
    @show typeof(G2), typeof(G4)
    println("Layer 1/$L")
    for i=1:L-1
        reg |> put(L, (i,i+1)=>G4)
    end
    for j=2:L
        println("Layer $j/$L")
        for i=1:L
            reg |> put(L, i=>G2)
        end
        for i=1:L-1
            reg |> put(L, (i,i+1)=>G4)
        end
    end
    sum(state(reg))
end

spinglass_yao(L::Int, jtype::Val; usecuda=false) = spinglass_yao(Float64, L, jtype; usecuda=usecuda)
function spinglass_yao(::Type{T}, L::Int, jtype::Val; usecuda=false) where T
    # Yao gates
    if usecuda
    	reg = ArrayReg(Tropical.(CuArrays.zeros(T, 1<<L)))
    else
    	reg = ArrayReg(ones(Tropical{T}, 1<<L))
    end
    _spinglass_yao(reg, L::Int, jtype::Val)
end
