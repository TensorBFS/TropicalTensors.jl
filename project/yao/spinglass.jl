using Test
using Yao
using CuYao
using LinearAlgebra
using CuArrays
CuArrays.allowscalar(false)
using TropicalTensors
using LuxurySparse
using DelimitedFiles

_get_J(::Val{:ferro}) = 1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()

function generate_J(type::Val, size...)
    J = zeros(size...)
    for i=1:length(J)
        J[i] = _get_J(type)
    end
    return J
end

G2(::Type{T}, J) where T = matblock(spinglass_bond_tensor(T(1.0)) |> LuxurySparse.staticize)
G4(::Type{T}, J) where T = matblock(Diagonal(spinglass_g4_tensor(T(1.0))) |> LuxurySparse.staticize)

function _spinglass_yao(reg::ArrayReg{B,Tropical{T}}, L::Int, J::AbstractVector) where {B,T}
    println("Layer 1/$L")
    k = 0
    for i=1:L-1
        k += 1
        reg |> put(L, (i,i+1)=>G4(T, J[k]))
    end
    for j=2:L
        println("Layer $j/$L")
        for i=1:L
            k += 1
            reg |> put(L, i=>G2(T, J[k]))
        end
        for i=1:L-1
            k += 1
            reg |> put(L, (i,i+1)=>G4(T, J[k]))
        end
    end
    sum(state(reg))
end

spinglass_yao(L::Int, J::AbstractVector; usecuda=false) = spinglass_yao(Float64, L, J; usecuda=usecuda)
function spinglass_yao(::Type{T}, L::Int, J::AbstractVector; usecuda=false) where T
    # Yao gates
    if usecuda
        reg = ArrayReg(CuArrays.ones(Tropical{T}, 1<<L))
    else
    	reg = ArrayReg(ones(Tropical{T}, 1<<L))
    end
    _spinglass_yao(reg, L, J)
end
