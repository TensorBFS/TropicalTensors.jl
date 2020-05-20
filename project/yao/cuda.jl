using CUDAnative
device!(1)
using Yao, CuYao
include("spinglass.jl")
include("datadump.jl")

L = 30
Js=load_J(L, Val(:randn))
@show Js
@show spinglass_yao(Float16, 30, Float16.(Js); usecuda=true)
@show spinglass_yao(Float32, 30, Float32.(Js); usecuda=true)
@show spinglass_yao(Float64, 30, Float64.(Js); usecuda=true)
