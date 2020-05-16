using TropicalTensors
using Viznet
using Yao
using TropicalYao

L = 18
Js = rand([-1, 1], 2L*(L-1))
hs = ones(Int, L^2)
@time solve(Spinglass(SquareLattice(L, L), Js, hs))

L = 22
Js = rand([-1, 1], 2L*(L-1))
hs = ones(Int, L^2)
@time gs = opt_config(Spinglass(SquareLattice(L, L), Js, hs))
