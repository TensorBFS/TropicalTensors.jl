using TropicalTensors
using DelimitedFiles
using Suppressor
using CuYao, CuArrays, CUDAnative
CUDAnative.device!(parse(Int, ARGS[1]))

using Random
Random.seed!(2)

function random_degeneracy(::Type{T}, lt; usecuda) where T
    sg = Spinglass(lt, rand(T[-1,1],length(sgbonds(lt))), zeros(T, length(lt)))
    res = solve(CountingTropical{T}, sg; usecuda=usecuda)
    return res.n, res.c
end

function main(::Type{T}, ::Type{LT}, Llist, nrepeat::Int; usecuda=false) where {T, LT}
    @suppress_err begin
        for (i, L) in enumerate(Llist)
            lt = LT(L, L)
            out = zeros(T, nrepeat, 2)
            for j = 1:nrepeat
                println("L=$L, $j")
                @time out[j,:] .= random_degeneracy(T, lt; usecuda=usecuda)
            end
            writedlm(joinpath(@__DIR__, "data", "$(name(LT))_degeneracy_L$L.dat"), out)
        end
    end
end

name(::Type{<:ChimeraLattice}) = "chimera"
name(::Type{<:SquareLattice}) = "square"
name(::Type{<:Cylinder}) = "cylinder"

if abspath(PROGRAM_FILE) == @__FILE__
    main(Int32, ChimeraLattice, 2:4, 1000; usecuda=false)
    main(Int32, ChimeraLattice, 5, 1000; usecuda=true)
    #main(Int32, SquareLattice, 8:4:16, 1000; usecuda=false)
    #main(Float64, SquareLattice, parse(Int, ARGS[2]), 1000; usecuda=true)
    #main(Float64, Cylinder, 4:4:12, 1000; usecuda=false)
    #main(Float64, Cylinder, [16, 20, 24], 1000; usecuda=true)
    #main(Int64, SquareLattice, [12 14], 1000; usecuda=false)
end
