using Test
using Yao
using CuYao
using LinearAlgebra
using CuArrays
CuArrays.allowscalar(false)
using TropicalTensors
using LuxurySparse
G2(::Type{T}, J) where T = matblock(spinglass_bond_tensor(T(J)) |> LuxurySparse.staticize)
G4(::Type{T}, J) where T = matblock(Diagonal(spinglass_g4_tensor(T(J))) |> LuxurySparse.staticize)
G16(::Type{T}, Js) where T = matblock(spinglass_g16_tensor(T.(Js)) |> LuxurySparse.staticize)

function red_reg(::Type{T}, Ly::Int, Js; usecuda=false) where T
    reg = _init_reg(T, Ly*4, usecuda)
    k = 0
    for i=1:Ly-1
        for j = 1:4
            reg |> put(4*Ly, (4*(i-1)+j,4*i+j)=>G4(T, Js[k+1]))
            k += 1
        end
    end
    for i=1:Ly
        reg |> put(4*Ly, (4i-3,4i-2,4i-1,4i)=>G16(T, Js[k+1:k+16]))
        k += 16
    end
    @assert k==length(Js)
    return reg
end

function chimera_yao(::Type{T}, Lx::Int, Ly::Int, J::AbstractVector; usecuda::Bool) where {T}
    println("Running Chimera, Lx = $Lx, Ly = $Ly, eltype = $T, usecuda = $usecuda.")
    nj_red = (Ly-1)*4 + 16 * Ly
    println("Layer 1/$Lx")
    k = 0
    reg = red_reg(T, Ly, J[k+1:k+nj_red]; usecuda=usecuda)
    k += nj_red

    for j=2:Lx
        println("Layer $j/$Lx")
        # BLACK
        for i=1:Ly*4
            reg |> put(Ly*4, i=>G2(T, J[k+1]))
            k += 1
        end

        # Contract with RED
        rr = red_reg(T, Ly, J[k+1:k+nj_red]; usecuda=usecuda)
        k += nj_red
        reg.state .*= rr.state
    end
    if length(J) !== k
        @warn "length of parameters is $(length(J)), used $k"
    end
    sum(state(reg))
end

chimera_yao(Lx::Int, Ly::Int, J::AbstractVector; usecuda=false) = chimera_yao(Float64, Lx, Ly, J; usecuda=usecuda)

function _init_reg(::Type{T}, L::Int, usecuda::Bool) where T
    if usecuda
        reg = ArrayReg(CuArrays.ones(Tropical{T}, 1<<L))
    else
    	reg = ArrayReg(ones(Tropical{T}, 1<<L))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Test
    @testset "test red reg" begin
        reg = red_reg(Float32, 3, randn(8+16*3); usecuda=false)
        @test reg isa ArrayReg{1,Tropical{Float32},Matrix{Tropical{Float32}}}
    end

    @testset "test Chimera" begin
        res = chimera_yao(Float32, 3, 3, ones(12*4 + 9*16); usecuda=false)
        @test res.n == 12*4 + 9*16
        res = chimera_yao(Float32, 3, 3, ones(12*4 + 9*16); usecuda=true)
        @test res.n == 12*4 + 9*16
    end
end
