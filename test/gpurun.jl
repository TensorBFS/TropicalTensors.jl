include("ttn2d.jl")

L = 30
e1 = @time (CuArrays.@sync gen_mine_2d(L, Val(:ferro), CC{:GPU,4}()))
@show e1
