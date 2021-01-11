using TropicalTensors
using TropicalTensors: getbonds, getvertices, rand_Js, rand_hs

function run(L)
    lt = CubicLattice(L, L, L)
    sg = Spinglass(lt, rand_Js(Float16, Randpm(), lt), rand_hs(Float16, Zero(), lt))
    res = @time solve(Tropical{Float16}, sg; usecuda=false)
end

const L = parse(Int, ARGS[1])
res = run(L)
println("The maximum energy is $(res.n).")
