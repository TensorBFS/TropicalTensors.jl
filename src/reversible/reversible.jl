module Reversible
using NiLang, NiLang.AD
using ..TropicalTensors
using LinearAlgebra

include("lib.jl")
include("ieinsum.jl")
include("blocks.jl")
include("instructs.jl")
include("tropicalgate.jl")

end
