using Test
using Yao
using LinearAlgebra
using TropicalTensors, TropicalTensors.Reversible
using LuxurySparse
using DelimitedFiles
using NiLang, NiLang.AD

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
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    @test spinglass_yao(Float32(0.0), reg, L, ones(Float32,180))[1] ≈ 180.0
    empty!(NiLang.GLOBAL_STACK)
end

function benchmarker(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    spinglass_yao(Float32(0.0), reg, L, ones(Float32, L*(L-1)*2))
    empty!(NiLang.GLOBAL_STACK)
end

function benchmarker2(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    _spinglass_yao(reg, L, ones(Float32, L*(L-1)*2))
end

include("../datadump.jl")
function opt_config(::Type{T}, L) where T
    Js = T.(load_J(L, Val(:randn)))
    reg = ArrayReg(ones(Tropical{T}, 1<<L))
    res = spinglass_yao(T(0.0), reg, L, Js)
    eng = res[1]
    gres = (~spinglass_yao)(GVar(eng, T(1.0)), GVar.(res[2:end])...)
    empty!(NiLang.GLOBAL_STACK)
    return eng, grad.(gres[end])
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

function assign_grid(L, g::AbstractVector) where {B,T<:Tropical}
    grid = zeros(Int, L, L)
    grid[1,1] = 1
    println("Layer 1/$L")
    k = 0
    for i=1:L-1
        k += identity(1)
        assign_one!(grid, (i,1), (i+1,1), g[k])
    end
    for j=2:L
        println("Layer $j/$L")
        for i=1:L
            k += identity(1)
            assign_one!(grid, (i,j-1), (i,j), g[k])
        end
        for i=1:L-1
            k += identity(1)
            assign_one!(grid, (i,j), (i+1,j), g[k])
        end
    end
    return grid
end

function assign_one!(grid, x, y, g)
    if grid[x...] == 0 && grid[y...] == 0
        error("")
    elseif grid[x...] == 0
        grid[x...] = sign(g)*grid[y...]
    elseif grid[y...] == 0
        grid[y...] = sign(g)*grid[x...]
    else
        @assert grid[y...] == sign(g)*grid[x...]
    end
end
