module TropicalTensors

using OMEinsum
using Requires

using TropicalNumbers, NiLog

export Tropical, TropicalF16, TropicalF32, TropicalF64
include("tensorlib.jl")
#include("tropical_mm.jl")
include("mislib.jl")

export Reversible
include("reversible/reversible.jl")

function __init__()
    @require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" include("cuda.jl")
end


end # module
