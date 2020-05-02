using CUDAnative
device!(7)
using Yao, CuYao
include("spinglass.jl")

@show spinglass_yao(Float32, 10, Val(:ferro); usecuda=true)
@show spinglass_yao(Float32, 10, Val(:ferro); usecuda=false)
