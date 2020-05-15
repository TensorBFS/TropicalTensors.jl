using TropicalTensors
using Test

@testset "dumpload" begin
    include("dumpload.jl")
end

@testset "yao" begin
    include("yao/yao.jl")
end

@testset "cuda" begin
    include("cuda.jl")
end
