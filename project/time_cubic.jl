using CUDA, CuYao
device!(parse(Int, ARGS[1]))
using TropicalTensors
using TropicalTensors: getbonds, getvertices, rand_Js, rand_hs

using BenchmarkTools
using Suppressor

function run()
    lt = CubicLattice(5, 5, 5)
    sg = Spinglass(lt, rand_Js(Float32, Randpm(), lt), rand_hs(Float32, Zero(), lt))
    @btime @suppress solve(Tropical{Float32}, $sg; usecuda=true)
    @btime @suppress solve(CountingTropical{Float32}, $sg; usecuda=true)
end


run()

