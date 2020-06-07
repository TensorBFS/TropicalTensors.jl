using TropicalTensors

lt = ChimeraLattice(8, 8)
sg = Spinglass(lt, randn(Int16, length(sgbonds(lt))), zeros(Int16, length(lt)))
res = @time solve(CountingTropical{T}, sg; usecuda=true)
res = @time solve(CountingTropical{T}, sg; usecuda=true)
@show res
