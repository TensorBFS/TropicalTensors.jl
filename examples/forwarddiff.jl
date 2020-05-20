using ForwardDiff, TropicalTensors
using CuYao

L = 28
Js = rand([-1, 1], 2L*(L-1))
@time solve(SquareLattice(L, L), Js, zeros(Int,L^2); usecuda=true)

L = 28
Js = rand([-1, 1], 2L*(L-1))
@time gs = ForwardDiff.gradient(x->solve(SquareLattice(L, L), convert.(eltype(x), Js),
    x; usecuda=true).n, zeros(Int, L^2))
