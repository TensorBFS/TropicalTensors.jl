using ForwardDiff, TropicalTensors
using CuYao

L = 5
Js = rand([-1, 1], 2L*(L-1))
gs = ForwardDiff.gradient(x->solve(SquareLattice(L, L), x,
    zeros(eltype(x), L^2); usecuda=true).n, Js)

using Compose
set_default_graphic_size(7cm, 7cm)
vizgrad_J(Spinglass(SquareLattice(L, L), Js, zeros(Int,L^2)), gs)
