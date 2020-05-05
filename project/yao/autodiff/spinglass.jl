using Test
using Yao
using LinearAlgebra
using TropicalTensors, TropicalTensors.Reversible
using LuxurySparse
using DelimitedFiles
using NiLang

@i function spinglass_yao(out!, reg::ArrayReg{B,T}, L::Int, J::AbstractVector) where {B,T<:Tropical}
    @safe println("Layer 1/$L")
    k ← 0
    for i=1:L-1
        k += identity(1)
        #reg |> put(L, (i,i+1)=>G4(T, J[k]))
        apply_G4!(reg, (i, i+1), J[k])
    end
    for j=2:L
        @safe println("Layer $j/$L")
        for i=1:L
            k += identity(1)
            #reg |> put(L, i=>G2(T, J[k]))
            apply_G2!(reg, i, J[k])
        end
        for i=1:L-1
            k += identity(1)
            #reg |> put(L, (i,i+1)=>G4(T, J[k]))
            apply_G4!(reg, (i, i+1), J[k])
        end
    end
    summed ← one(T)
    isum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(J)
    summed → one(T)
end

@testset "yao" begin
    L = 10
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    @test check_inv(spinglass_yao, (Float32(0.0), reg, L, ones(Float32, 180)))
    @test spinglass_yao(Float32(0.0), reg, L, ones(Float32,180))[1] ≈ 180.0
    empty!(NiLang.GLOBAL_STACK)
end

function benchmarker(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    spinglass_yao(Float32(0.0), reg, L, ones(L*(L-1)*2))
end

function benchmarker2(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    _spinglass_yao(reg, L, ones(L*(L-1)*2))
end
