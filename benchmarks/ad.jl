@assert ARGS[1] in ["chimera", "square"]
@assert ARGS[2] in ["forwarddiff", "nilang"]

const AD = ARGS[2]

using CUDAnative, CuArrays
device!(parse(Int, get(ARGS, 3, "0")))
using CuYao
using DelimitedFiles, TropicalTensors
using BenchmarkTools
using ForwardDiff

if ARGS[1] == "chimera"
    const LATTICE = ChimeraLattice
else
    const LATTICE = SquareLattice
end

suite = BenchmarkGroup()

function bfunc(::Type{T}, L::Int) where T
    lt = LATTICE(L, L)
    if T <: Integer
        Js = rand(T[-1, 1], length(sgbonds(lt)))
    else
        Js = randn(T, length(sgbonds(lt)))
    end
    hs = zeros(T, length(lt))
    sg = Spinglass(lt, Js, hs)
    if AD == "forwarddiff"
        CuArrays.@sync ForwardDiff.gradient(x->solve(SquareLattice(L, L),
            convert.(eltype(x), Js),
            x; usecuda=true).n, hs)
    else
        opt_config(sg)
    end
end

if LATTICE === SquareLattice
    const RANGE = 4:2:25
    const TYPE = Float32
else
    const RANGE = 2:6
    const TYPE = Int16
end

for L = RANGE
    suite[L] = @benchmarkable bfunc($TYPE, $L)
end

println("loaded")

tune!(suite)
res = run(suite)

function analyze_res(res)
    times = zeros(length(res))
    for (k,L) = enumerate(RANGE)
        times[k] = minimum(res[L].times)
    end
    return times
end

times = analyze_res(res)
fname = joinpath(@__DIR__, "data", "bench_$(LATTICE)_$AD.dat")
println("Writing benchmark results to file: $fname.")
writedlm(fname, times)
