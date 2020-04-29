using Yao, CuYao
include("spinglass.jl")

@time spinglass_yao(Float32, 32, Val(:ferro); usecuda=true)
