module Reversible
using NiLang, NiLang.AD
using ..TropicalTensors
using LinearAlgebra
using TropicalNumbers
using NiLog

include("blocks.jl")
include("instructs.jl")
include("tropicalgate.jl")

end
