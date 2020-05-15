using DelimitedFiles
using Random

export dump_params, load_params

# chimera
function dump_params(folder, sg::Spinglass{LT}, jt::JT, ht::HT, seed::Int) where {LT,JT<:Jtype, HT<:Jtype}
    jfilename = joinpath(folder, _gen_jfname(sg, jt, seed))
    hfilename = joinpath(folder, _gen_hfname(sg, ht, seed))
    open(jfilename, "w") do f
        write(f, "# $JT\n")
        for ((i,j), J) in zip(sgbonds(sg.lattice), sg.Js)
            write(f, "$i $j $J\n")
        end
    end
    open(hfilename, "w") do f
        write(f, "# $HT\n")
        for (i, h) in zip(vertices(sg.lattice), sg.hs)
            write(f, "$i $h\n")
        end
    end
    return sg
end

function _gen_jfname(sg::Spinglass{LT}, ::JT, seed::Int) where {LT,JT, HT}
    sz = size(sg.lattice)
    return "$(LT)_J$(JT)_Lx$(sz[end-1])_Ly$(sz[end])_seed$(seed).dat"
end

function _gen_hfname(sg::Spinglass{LT}, ::HT, seed::Int) where {LT,JT, HT}
    sz = size(sg.lattice)
    return "$(LT)_h$(HT)_Lx$(sz[end-1])_Ly$(sz[end])_seed$(seed).dat"
end

function load_params(folder, sg::Spinglass{LT,T}, jt::JT, ht::HT, seed::Int) where {LT,T,JT<:Jtype, HT<:Jtype}
    jfilename = joinpath(folder, _gen_jfname(sg, jt, seed))
    hfilename = joinpath(folder, _gen_hfname(sg, ht, seed))
    open(jfilename, "r") do f
        readline(f)
        for (k, (i,j)) in enumerate(sgbonds(sg.lattice))
            line = strip(readline(f))
            sg.Js[k] = parse(T, split(line, ' ')[3])
        end
    end
    open(hfilename, "r") do f
        readline(f)
        for i in vertices(sg.lattice)
            line = strip(readline(f))
            sg.hs[i] = parse(T, split(line, ' ')[2])
        end
    end
    return sg
end
