using CUDAnative
device!(2)
include("ttn2d.jl")

L = 30
e1 = @time (CuArrays.@sync gen_mine_2d(L, Val(:ferro), CC{:CPU,1}()))
@show e1
