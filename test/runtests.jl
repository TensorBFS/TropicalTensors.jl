using TropicalTensors
using Test

@testset "dumpload" begin
    include("dumpload.jl")
end

@testset "yao" begin
    include("yao.jl")
end

@testset "reversible" begin
    include("reversible/reversible.jl")
end

@testset "cuda" begin
    if Base.find_package("CuYao") != nothing
        include("cuda.jl")
    end
end
