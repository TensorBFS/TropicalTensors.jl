using CUDAnative
device!(2)
include("ttn2d.jl")

L = 28
println("Running for L = $L")
e1 = @time (CuArrays.@sync gen_mine_2d(L, Val(:ferro), CC{:GPU,1}()))
@show e1
