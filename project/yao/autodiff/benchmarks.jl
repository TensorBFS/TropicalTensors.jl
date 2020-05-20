using BenchmarkTools
include("lib.jl")

N = 50
x = Tropical.(randn(N, N))
y = Tropical.(randn(N, N))
out2 = Tropical.(zeros(N,N))

@benchmark $x * $y
@benchmark igemm!($out2, $x, $y)


N = 500
x = Tropical.(randn(N, N))
v = Tropical.(randn(N))
out2 = Tropical.(zeros(N))

using StaticArrays
m = MMatrix{16,16}(Tropical.(randn(16,16)))
v = Tropical.(randn(16))
@benchmark $m * $v
bk = zeros(Bool, length(v))
out = zero(v)
@benchmark @instr igemv!($out, $m, $v, $bk)

#### instructs  ###########
using NiLang, NiLang.AD
using TropicalTensors
using LinearAlgebra, Yao
using StaticArrays

include("lib.jl")
include("instructs.jl")

using BenchmarkTools
g4 = Diagonal(MMatrix{4,4}(spinglass_g4_tensor(1.5)))
reg = ArrayReg(randn(1<<18) .|> Tropical)
@benchmark i_instruct!($(copy(vec(reg.state))), $g4, (3, 7), (), ())
@benchmark i_instruct!($(copy(GVar.(vec(reg.state)))), $(GVar(g4)), (3, 7), (), ())
@benchmark $(copy(reg)) |> $(put(18, (3, 7)=>matblock(g4)))

g4 = MMatrix{4,4}(Tropical.(randn(4,4)))
@benchmark (i_instruct!($(copy(vec(reg.state))), $g4, (3, 7), (), ()); empty!(NiLang.GLOBAL_STACK))
@benchmark (i_instruct!($(copy(GVar.(vec(reg.state)))), $(GVar(g4)), (3, 7), (), ()); empty!(NiLang.GLOBAL_STACK))
@benchmark copy($reg) |> $(put(18, (3, 7)=>matblock(g4)))


################ naive einsum #############
using BenchmarkTools

N = 50
a = Tropical.(randn(N, N))
b = Tropical.(randn(N, N))
c = Tropical.(zeros(N, N))
y = Tropical.(zeros(N, N, N))
@benchmark naive_einsum!(ein"ab,bc->ac", ($a,$b), c) seconds = 1
@benchmark naive_einsum!(((1,2), (1,3), (1,4)), ($a,$a,$a), (2,3,4), $y) seconds = 1
