using TropicalTensors
using Test
using OMEinsum

@testset "contract" begin
    include("mislib.jl")
end

@testset "tensorlib" begin
    include("tensorlib.jl")
end

@testset "reversible" begin
    include("reversible/reversible.jl")
end
