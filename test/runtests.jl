using TropicalTensors
using Test

@testset "dumpload" begin
    include("dumpload.jl")
end

@testset "yao" begin
    include("yao/yao.jl")
end

@testset "yao reversible" begin
    include("yao/reversible/square.jl")
    include("yao/reversible/chimera.jl")
end

@testset "cuda" begin
    include("cuda.jl")
end
