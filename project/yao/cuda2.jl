using CUDAnative
device!(7)
using Yao, CuYao
include("chimera.jl")
include("datadump.jl")

L = 6
@time chimera_yao(Float32, L, L, load_JCs(L, L, Val(:randn)); usecuda=true)
#@time chimera_yao(Float32, L, L, load_JCs(L, L, Val(:randn)); usecuda=true)
