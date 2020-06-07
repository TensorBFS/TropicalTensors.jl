using TropicalTensors, DelimitedFiles

function run_chimera44(::Type{T}; nrepeat, usecuda=false) where T
    holes = [9, 11, 12, 20, 24, 28, 29, 30,
        32, 59,
        70, 89, 90, 91,
        99, 111, 118, 121, 123, 127] .+ 1
    lt = ChimeraLattice(4, 4)
    bonds = sgbonds(lt)
    @assert length(holes) == 20
    out = zeros(T, 2, nrepeat)

    for j=1:nrepeat
        # set some of `J`s to 0
        Js = rand(T[-1, 1], length(bonds))
        hs = zeros(T, 16*8)
        for i=1:length(Js)
            if (bonds[i][1] in holes) || (bonds[i][2] in holes)
                Js[i] = 0
            end
        end
        sg = Spinglass(lt, Js, hs)
        #display(sg)
        res = solve(CountingTropical{T}, sg; usecuda=usecuda)
        out[:,j] .= (res.n, res.c)
        ntimes = res.c % (2^(length(holes)+1))
        if ntimes != 0
            error("got ($(res.c)% $(2^(length(holes)+1)))")
        end
    end
    writedlm(joinpath(@__DIR__, "data", "chimera44_degeneracy.dat"), out)
    return out
end

run_chimera44(Float32; nrepeat=1000, usecuda=false)
