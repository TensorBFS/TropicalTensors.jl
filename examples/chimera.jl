using ForwardDiff, TropicalTensors
using CuYao

L = 8
lt = ChimeraLattice(L, L)
Js = rand([Int16(-1), Int16(1)], length(sgbonds(lt)))
hs = zeros(Int16, length(lt))
@time solve(lt, Js, hs; usecuda=true)
