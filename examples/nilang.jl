using TropicalTensors
using Viznet
using Yao
using TropicalYao

L = 18
Js = rand([-1, 1], 2L*(L-1))
hs = ones(Int, L^2)
@time solve(Spinglass(SquareLattice(L, L), Js, hs))

L = 7
lt = ChimeraLattice(L, L)
Js = rand(Int16[-1, 1], length(sgbonds(lt)))
hs = ones(Int16, length(lt))
tofile("chimera77.dat", @time gs = opt_config(Spinglass(lt, Js, hs)))
