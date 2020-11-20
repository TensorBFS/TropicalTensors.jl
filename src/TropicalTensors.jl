module TropicalTensors

using Viznet
using Requires

using TropicalNumbers

export Tropical, CountingTropical
export Reversible

project_relative_path(xs...) = normpath(joinpath(dirname(dirname(pathof(@__MODULE__))), xs...))

include("Core.jl")
include("dumpload.jl")
include("yao.jl")
include("tensorcontract.jl")
#export Reversible
include("reversible/reversible.jl")
include("viz.jl")

function __init__()
    @require CuYao = "b48ca7a8-dd42-11e8-2b8e-1b7706800275" include("cuda.jl")
    @require Pluto = "c3e4b0f8-55cb-11ea-2926-15256bba5781" include("pluto.jl")
end

end # module
