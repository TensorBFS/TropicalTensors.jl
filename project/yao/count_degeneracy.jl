using TropicalTensors
using DelimitedFiles
using CuYao, CuArrays, CUDAnative
CUDAnative.device!(1)

function random_chimera_degeneracy(Lx::Int, Ly::Int; usecuda)
    lt = ChimeraLattice(Lx, Ly)
    sg = Spinglass(lt, rand(Int32[-1,1],length(sgbonds(lt))), zeros(Int32, length(lt)))
    res = solve(CountingTropical{Int32}, sg; usecuda=usecuda)
    return res.c
end

function main(Llist, nrepeat::Int; usecuda=false)
    for (i, L) in enumerate(Llist)
        out = zeros(Int, nrepeat)
        for j = 1:nrepeat
            println("L=$L, $j")
            @time out[j] = random_chimera_degeneracy(L, L; usecuda=usecuda)
        end
        writedlm(joinpath(@__DIR__, "data", "chimera_degeneracy_L$L.dat"), out)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(6, 1000; usecuda=true)
end
