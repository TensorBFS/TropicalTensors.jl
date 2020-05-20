using TropicalTensors
using Test

@testset "dumpload" begin
    include("dumpload.jl")
end

@testset "yao" begin
    include("yao/yao.jl")
end

@testset "reversible square" begin
    include("yao/reversible/square.jl")
end

@testset "reversible chimera" begin
    include("yao/reversible/chimera.jl")
end

@testset "cuda" begin
    include("cuda.jl")
end
