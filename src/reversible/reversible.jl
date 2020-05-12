module Reversible
using NiLang, NiLang.AD
using ..TropicalTensors
using LinearAlgebra
using TropicalNumbers

include("lib.jl")
include("ieinsum.jl")
include("blocks.jl")
include("instructs.jl")
include("tropicalgate.jl")

end
