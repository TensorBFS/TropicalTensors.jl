using DelimitedFiles
using CUDA, CuYao
device!(parse(Int, ARGS[1]))

using Random
using TropicalTensors

function run(::Type{T}, n::Int; nrepeat, usecuda) where T
    saveto = joinpath(@__DIR__, "potts_n$(n)_elsl.dat")
    elsl = zeros(T, 3, nrepeat)
    println("number of repeatition is $nrepeat.")
    for i=1:nrepeat
        lt = SquareLattice(n, n)
        t = @elapsed res = solve_potts(CountingTropical{T}, Val(3), lt, TropicalTensors.potts_randpm_J(lt); usecuda=usecuda)
        println("Δt = $t, maximum energy = $(res.n), degeneracy = $(res.c)")
        elsl[1,i] = res.n/n
        elsl[2,i] = log(res.c)/n
        elsl[3,i] = t
    end
    writedlm(saveto, elsl)
    return elsl
end

L = parse(Int, ARGS[2])
run(Float32, L; nrepeat=100, usecuda=true)