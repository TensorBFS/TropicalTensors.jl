using CUDAnative
device!(2)
using Yao, CuYao
include("chimera.jl")
include("datadump.jl")

L = 8
#@time chimera_yao(Float32, L, L, load_JC(L, L, Val(:randn)); usecuda=false)
res = @time chimera_yao(Float16, L, L, load_JC(L, L, Val(:ferro)); usecuda=true)
@show res
#@time chimera_yao(Float32, L, L, load_JCs(L, L, Val(:randn)); usecuda=true)
#@time chimera_yao(Float32, L, L, load_JCs(L, L, Val(:randn)); usecuda=true)
