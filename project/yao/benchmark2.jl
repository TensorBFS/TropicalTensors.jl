using CUDAnative
using DelimitedFiles
device!(2)

include("chimera.jl")
include("datadump.jl")

# dump_Js(Val(:randn))
#using Test
#J = load_J(10, Val(:ferro))
#@test spinglass_yao(10, J; usecuda=false).n == 180
#@test spinglass_yao(10, J; usecuda=true).n == 180

using BenchmarkTools
suite = BenchmarkGroup()
suite["GPU"] = BenchmarkGroup()
for L = 3:8
    suite["GPU"][L] = @benchmarkable (CuArrays.@sync chimera_yao(Float16, $L, $L, $(load_JC(L, L, Val(:randn))); usecuda=true))
    #suite["CPU"][L] = @benchmarkable (CuArrays.@sync chimera_yao(Float16, $L, $(load_J(L, Val(:randn))); usecuda=false))
end

println("loaded")

tune!(suite)
res = run(suite)

function analyze_res(res)
    times = zeros(length(res["GPU"]))
    for (k,L) = enumerate(3:8)
        times[k] = minimum(res["GPU"][L].times)
    end
    return times
end

times = analyze_res(res)
fname = joinpath(@__DIR__, "data", "bench_chimera_gpu.dat")
println("Writing benchmark results to file: $fname.")
writedlm(fname, times)
