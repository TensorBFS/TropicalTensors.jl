module TropicalTensors

using OMEinsum
using Requires

include("tropical_number.jl")
include("tensorlib.jl")
#include("tropical_mm.jl")
include("mislib.jl")

export Reversible
include("reversible/reversible.jl")

function __init__()
    @require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" include("cuda.jl")
end


end # module
