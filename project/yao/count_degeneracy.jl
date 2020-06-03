using TropicalTensors
using DelimitedFiles

function random_chimera_degeneracy(Lx::Int, Ly::Int)
    lt = ChimeraLattice(Lx, Ly)
    sg = Spinglass(lt, rand(Int32[-1,1], 12*4 + 9*16), zeros(Int32, 72))
    res = solve(CountingTropical{Int32}, sg; usecuda=false)
    return res.c
end

function main(Llist, nrepeat::Int)
    out = zeros(Int, nrepeat, length(Llist))
    for (i, L) in enumerate(Llist)
        for j = 1:nrepeat
            out[j, i] = random_chimera_degeneracy(L, L)
        end
    end
    writedlm(joinpath(@__DIR__, "data", "chimera_degeneracy.dat"), out)
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(2:6, 1000)
end
