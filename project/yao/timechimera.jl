using TropicalTensors
using CuYao, CUDAnative, CuArrays
CUDAnative.device!(6)
using Suppressor

const T = Float16
lt = ChimeraLattice(8, 8)
sg = Spinglass(lt, rand(T[1], length(sgbonds(lt))), zeros(T, length(lt)))
res = @time begin
    @suppress solve(Tropical{T}, sg; usecuda=true)
end
res = @time begin
    @suppress solve(Tropical{T}, sg; usecuda=true)
end
@show res
