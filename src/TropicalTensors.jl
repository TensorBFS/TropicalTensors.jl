module TropicalTensors

using Viznet
using Requires

using TropicalNumbers

export Tropical, CountingTropical
export Reversible

include("Core.jl")
include("dumpload.jl")
include("yao.jl")
include("tensorcontract.jl")
#export Reversible
include("reversible/reversible.jl")
include("viz.jl")
include("pluto.jl")

function __init__()
    @require CuYao = "b48ca7a8-dd42-11e8-2b8e-1b7706800275" include("cuda.jl")
end

end # module
