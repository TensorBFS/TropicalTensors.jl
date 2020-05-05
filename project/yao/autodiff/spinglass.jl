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
    @test spinglass_yao(0.0, reg, L, ones(180))[1] ≈ 180.0
end
