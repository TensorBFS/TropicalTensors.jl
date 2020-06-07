using TropicalTensors, DelimitedFiles
using Suppressor

const holes = [9, 11, 12, 20, 24, 28, 29, 30,
    32, 59,
    70, 89, 90, 91,
    99, 111, 118, 121, 123, 127] .+ 1

function load_instance(filename::String)
    data = readdlm(filename, comments=true, comment_char='#')
    lt = ChimeraLattice(4, 4)
    vertices = sgvertices(lt)
    map = setdiff(1:128, holes)
    res = Dict{Tuple{Int, Int}, Int}()
    for i = 1:size(data, 1)
        a, b, J = data[i,:]
        a = Int(a)
        b = Int(b)
        res[map[a], map[b]] = Int(J)
    end
    res
end

data = load_instance(joinpath(@__DIR__, "anc",
    "instances", "Benchmarking--108q--15us--02-Aug-2012-13-55-20.txt"))

minfirst(x) = x[1] > x[2] ? (x[2], x[1]) : x

function run_chimera44(::Type{T}; filenames, usecuda=false) where T
    lt = ChimeraLattice(4, 4)
    bonds = sgbonds(lt)
    @assert length(holes) == 20
    out = zeros(T, 2, length(filenames))

    for j=1:length(filenames)
        @show j, filenames[j]
        # set some of `J`s to 0
        Js = rand(T[-1, 1], length(bonds))
        hs = zeros(T, 16*8)
        d = load_instance(filenames[j])
        for i=1:length(Js)
            if (bonds[i][1] in holes) || (bonds[i][2] in holes)
                Js[i] = 0
            else
                Js[i] = d[minfirst(bonds[i])]
            end
        end
        sg = Spinglass(lt, Js, hs)
        #display(sg)
        res = @suppress begin
            solve(CountingTropical{T}, sg; usecuda=usecuda)
        end
        out[:,j] .= (res.n, res.c)
        ntimes = res.c % (2^(length(holes)+1))
        if ntimes != 0
            error("got ($(res.c)% $(2^(length(holes)+1)))")
        end
    end
    writedlm(joinpath(@__DIR__, "data", "chimera44_degeneracy.dat"), out)
    return out
end

filenames = readdir(joinpath(@__DIR__, "anc", "instances"), join=true)
run_chimera44(Float32; usecuda=false, filenames=filenames)

function load_success(filename)
    res = open(filename) do f
        readlines(f)
    end
    map(x->parse(Float64, split(x, " ")[2]), res)
end

suc = load_success(joinpath(@__DIR__, "anc", "success.txt"))
writedlm(joinpath(@__DIR__, "data", "success.dat"), suc)
