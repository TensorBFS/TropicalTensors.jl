using CUDAnative
device!(7)
using Yao, CuYao
include("chimera.jl")

@time chimera_yao(Float32, 8, 8, ones(8*7*8 + 64*16); usecuda=true)
