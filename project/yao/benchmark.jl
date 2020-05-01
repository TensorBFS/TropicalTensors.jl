using CUDAnative
using DelimitedFiles
device!(6)

include("spinglass.jl")

function dump_Js(jtype::Val{JT}) where JT
    for i=4:2:32
        J = generate_J(jtype, i*(i-1)*2)
        writedlm(joinpath(@__DIR__, "data", "J_$(JT)_L$i.dat"), J)
    end
end

function load_J(L::Int, jtype::Val{JT}) where JT
    vec(readdlm(joinpath(@__DIR__, "data", "J_$(JT)_L$L.dat")))
end

# dump_Js(Val(:randn))
#using Test
#J = load_J(10, Val(:ferro))
#@test spinglass_yao(10, J; usecuda=false).n == 180
#@test spinglass_yao(10, J; usecuda=true).n == 180

using BenchmarkTools
suite = BenchmarkGroup()
suite["GPU"] = BenchmarkGroup()
for L = 4:2:10
    suite["GPU"][L] = @benchmarkable spinglass_yao($L, $(load_J(L, Val(:randn))); usecuda=true)
    #suite["CPU"][L] = @benchmarkable spinglass_yao(L, $(load_J(L, Val(:randn))); usecuda=false)
end

tune!(suite)
res = run(suite)

function analyze_res(res)
    times = zeros(length(res["GPU"]))
    for (k,L) = enumerate(4:2:10)
        times[k] = minimum(res["GPU"][L].times)
    end
    return times
end

times = analyze_res(res)
fname = joinpath(@__DIR__, "data", "bench_spingglass_gpu.dat")
println("Writing benchmark results to file: $fname.")
writedlm(fname, times)
