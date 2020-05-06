using Test
using Yao
using LinearAlgebra
using TropicalTensors, TropicalTensors.Reversible
using LuxurySparse
using DelimitedFiles
using NiLang, NiLang.AD

@i function spinglass_yao(out!, reg::ArrayReg{B,T}, L::Int, J::AbstractVector, A_STACK, B_STACK) where {B,T<:Tropical}
    @invcheckoff begin
    @safe println("Layer 1/$L, stack size: $(A_STACK.top) & $(B_STACK.top)")
    k ← 0
    for i=1:L-1
        k += identity(1)
        apply_G4!(reg, (i, i+1), J[k], A_STACK)
    end
    for j=2:L
        @safe println("Layer $j/$L, stack size: $(A_STACK.top) & $(B_STACK.top)")
        @routine begin
            for i=1:L
                k += identity(1)
                apply_G2!(reg, i, J[k], A_STACK)
            end
        end
        # bookup current state and loss
        incstack!(B_STACK)
        for i = 1:length(reg.state)
            @inbounds muleq(B_STACK[i], reg.state[i])
        end

        # clean up `NiLang.GLOBAL_STACK`
        ~@routine
        # but wait, `k` should not be uncomputed
        k += identity(L)

        # store the stored state
        for i=1:length(reg.state)
            @inbounds NiLang.SWAP(reg.state[i], B_STACK[i])
        end

        for i=1:L-1
            k += identity(1)
            apply_G4!(reg, (i, i+1), J[k], A_STACK)
        end
    end
    summed ← one(T)
    isum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(J)
    summed → one(T)
    end
end

@i function spinglass_yao_v0(out!, reg::ArrayReg{B,T}, L::Int, J::AbstractVector, REG_STACK) where {B,T<:Tropical}
    @invcheckoff begin
    @safe println("Layer 1/$L")
    k ← 0
    for i=1:L-1
        k += identity(1)
        apply_G4!(reg, (i, i+1), J[k], REG_STACK)
    end
    for j=2:L
        @safe println("Layer $j/$L")
        for i=1:L
            k += identity(1)
            apply_G2!(reg, i, J[k], REG_STACK)
        end
        for i=1:L-1
            k += identity(1)
            apply_G4!(reg, (i, i+1), J[k], REG_STACK)
        end
    end
    summed ← one(T)
    isum(summed, reg.state)
    NiLang.SWAP(summed.n, out!)
    k → length(J)
    summed → one(T)
    end
end

@testset "yao" begin
    L = 10
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    A = stack4reg(reg, L)
    B = stack4reg(reg, L-1)
    @test check_inv(spinglass_yao, (Float32(0.0), reg, L, ones(Float32, 180), A, B))
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    A = stack4reg(reg, L)
    B = stack4reg(reg, L-1)
    @test spinglass_yao(Float32(0.0), reg, L, ones(Float32,180), A, B)[1] ≈ 180.0
    empty!(NiLang.GLOBAL_STACK)
end

function benchmarker(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    A = stack4reg(reg, L)
    B = stack4reg(reg, L-1)
    spinglass_yao(Float32(0.0), reg, L, ones(Float32, L*(L-1)*2), A, B)
    println("cached array size: $(count(x->x isa AbstractArray, NiLang.GLOBAL_STACK))")
    empty!(NiLang.GLOBAL_STACK)
end

function benchmarker1(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    spinglass_yao_v0(Float32(0.0), reg, L, ones(Float32, L*(L-1)*2))
    println("cached array size: $(count(x->x isa AbstractArray, NiLang.GLOBAL_STACK))")
    empty!(NiLang.GLOBAL_STACK)
end

function benchmarker2(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    _spinglass_yao(reg, L, ones(Float32, L*(L-1)*2))
end

function benchmarker4(L)
    reg = ArrayReg(Tropical.(map(x->ForwardDiff.Dual(x, zero(x)), zeros(1<<L))))
    Js = ones(Float32, L*(L-1)*2)
    dJs = map(i->ForwardDiff.Dual(Js[i], i==1 ? one(Js[i]) : zero(Js[i])), 1:length(Js))
    _spinglass_yao(reg, L, dJs)
end

# fail
function benchmarker5(L)
    function f(Js)
        reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
        _spinglass_yao(reg, L, Js).n
    end
    Zygote.gradient(f, ones(Float32, L*(L-1)*2))
end

include("../datadump.jl")

function opt_config(::Type{T}, L; jtype=:randn) where T
    Js = T.(load_J(L, Val(jtype)))
    reg = ArrayReg(ones(Tropical{T}, 1<<L))
    A = stack4reg(reg, L)
    B = stack4reg(reg, L-1)
    eng, reg, L, Js, A, B = spinglass_yao(T(0.0), reg, L, Js, A, B)
    println("size of global stack:", length(NiLang.GLOBAL_STACK))
    gres = (~spinglass_yao)(GVar(eng, T(1)), GVar(reg), L, GVar.(Js, zero(Js)), GVar(A), GVar(B))
    empty!(NiLang.GLOBAL_STACK)
    return eng, grad.(gres[end-2])
end

G2(::Type{T}, J) where T = matblock(spinglass_bond_tensor(T(J)) |> LuxurySparse.staticize)
G4(::Type{T}, J) where T = matblock(Diagonal(spinglass_g4_tensor(T(J))) |> LuxurySparse.staticize)

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

function print_grid(grid)
    M, N = size(grid)
    for i=1:M
        for j=1:N
            dot = grid[i,j] == 1 ? "⚫" : "⚪"
            print(" $dot ")
        end
        println()
    end
end
