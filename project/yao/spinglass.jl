using Test
using Yao
using CuYao
using CuArrays
CuArrays.allowscalar(false)
using TropicalTensors

_get_J(::Val{:ferro}) = 1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()

function _spinglass_yao(reg, L::Int, jtype::Val)
    G2 = matblock(spinglass_bond_tensor(1.0))
    G4 = matblock(spinglass_g4_tensor(1.0))
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
    sum(statevec(reg))
end

spinglass_yao(L::Int, jtype::Val; usecuda=false) = spinglass_yao(L, jtype; usecuda=usecuda)
function spinglass_yao(::Type{T}, L::Int, jtype::Val; usecuda=false) where T
    # Yao gates
    reg = ArrayReg(ones(Tropical{T}, 1<<L))
    if usecuda
        reg = reg |> cu
    end
    _spinglass_yao(reg, L::Int, jtype::Val)
end
