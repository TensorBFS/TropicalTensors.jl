include("chimera.jl")
include("datadump.jl")

L = 4
chimera_yao(Float32, L, L, load_JC(L, L, Val(:randpm)); usecuda=false)
