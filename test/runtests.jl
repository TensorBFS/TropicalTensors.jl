using TropicalTensors
using Test

@testset "dumpload" begin
    include("dumpload.jl")
end

@testset "yao" begin
    include("yao.jl")
end

@testset "reversible square" begin
    include("reversible/square.jl")
end

@testset "reversible chimera" begin
    include("reversible/chimera.jl")
end

@testset "reversible second_neighbor" begin
    include("reversible/second_neighbor.jl")
end

@testset "cuda" begin
    include("cuda.jl")
end
