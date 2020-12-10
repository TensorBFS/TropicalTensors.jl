using Random
using TropicalTensors
using DelimitedFiles
using CuYao, CUDA
CUDA.device!(parse(Int, ARGS[1]))

function random_degeneracy(::Type{T}, lt; usecuda) where T
    sg = Spinglass(lt, rand(T[-1,1],length(sgbonds(lt))), zeros(T, length(lt)))
    res = solve(CountingTropical{T}, sg; usecuda=usecuda)
    return res.n, res.c
end

function main(::Type{T}, ::Type{LT}, Llist, nrepeat::Int; usecuda=false) where {T, LT}
    println("The number of repeatition is $nrepeat.")
    for (i, L) in enumerate(Llist)
        lt = LT(L, L)
        out = zeros(T, nrepeat, 2)
        for j = 1:nrepeat
            println("L=$L, $j")
            @time out[j,:] .= random_degeneracy(T, lt; usecuda=usecuda)
        end
        fname = joinpath(@__DIR__, "$(name(LT))_degeneracy_L$L.dat")
        writedlm(fname, out)
    end
end

name(::Type{<:ChimeraLattice}) = "chimera"
name(::Type{<:SquareLattice}) = "square"
name(::Type{<:Cylinder}) = "cylinder"

if abspath(PROGRAM_FILE) == @__FILE__
    Random.seed!(2)
    nrepeat = 100
    L = parse(Int, ARGS[2])
    main(Float32, Cylinder, [L], nrepeat; usecuda=true)
end
