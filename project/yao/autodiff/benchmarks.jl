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
