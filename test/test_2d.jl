include("ttn2d.jl")

L = 16
e1 = gen_mine_2d(L, Val(:ferro), CC{:GPU,4}())
e2 = gen_mine_2d(L, Val(:ferro), CC{:CPU,4}())
e3 = gen_mine_2d(L, Val(:ferro), CC{:GPU,1}())
e4 = gen_mine_2d(L, Val(:ferro), CC{:CPU,1}())
using Test
@test Array(e1) == e2 == Array(e3) == e4

using BenchmarkTools
@benchmark gen_mine_2d($L, $(Val(:ferro)), CC{:CPU,1}())

using CuArrays
e = gen_mine_2d(L, Val(:ferro), CC{:GPU,1}())
@benchmark (CuArrays.@sync gen_mine_2d($L, $(Val(:ferro)), CC{:GPU,1}()))
