using TropicalTensors
using Test
using OMEinsum
using TropicalYao

@testset "contract" begin
    include("mislib.jl")
end

@testset "dumpload" begin
    include("dumpload.jl")
end
