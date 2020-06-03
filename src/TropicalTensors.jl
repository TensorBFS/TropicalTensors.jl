module TropicalTensors

using Viznet
using Requires

using TropicalNumbers

export Tropical, CountingTropical

include("Core.jl")
include("dumpload.jl")
include("yao.jl")
#export Reversible
include("reversible/reversible.jl")

function __init__()
    @require CuYao = "b48ca7a8-dd42-11e8-2b8e-1b7706800275" include("cuda.jl")
end

end # module
