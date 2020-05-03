using CUDAnative
using DelimitedFiles
device!(6)

include("spinglass.jl")
include("datadump.jl")

# dump_Js(Val(:randn))
#using Test
#J = load_J(10, Val(:ferro))
#@test spinglass_yao(10, J; usecuda=false).n == 180
#@test spinglass_yao(10, J; usecuda=true).n == 180

using BenchmarkTools
suite = BenchmarkGroup()
suite["GPU"] = BenchmarkGroup()
for L = 4:2:32
    suite["GPU"][L] = @benchmarkable (CuArrays.@sync spinglass_yao(Float32, $L, $(load_J(L, Val(:randn))); usecuda=true))
    #suite["CPU"][L] = @benchmarkable (CuArrays.@sync spinglass_yao(Float32, $L, $(load_J(L, Val(:randn))); usecuda=false))
end

println("loaded")

tune!(suite)
res = run(suite)

function analyze_res(res)
    times = zeros(length(res["GPU"]))
    for (k,L) = enumerate(4:2:32)
        times[k] = minimum(res["GPU"][L].times)
    end
    return times
end

times = analyze_res(res)
fname = joinpath(@__DIR__, "data", "bench_spinglass_gpu.dat")
println("Writing benchmark results to file: $fname.")
writedlm(fname, times)
