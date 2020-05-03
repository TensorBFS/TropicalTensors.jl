using DelimitedFiles

_get_J(::Val{:ferro}) = 1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()

function generate_J(type::Val, size...)
    J = zeros(size...)
    for i=1:length(J)
        J[i] = _get_J(type)
    end
    return J
end

# chimera
function dump_JCs(jtype::Val{JT}) where JT
    for i=3:8
        J = generate_J(jtype, i*(i-1)*8 + i^2*16)
        writedlm(joinpath(@__DIR__, "data", "JC_$(JT)_Lx$(i)_Ly$(i).dat"), J)
    end
end

function load_JC(Lx::Int, Ly::Int, jtype::Val{JT}) where JT
    vec(readdlm(joinpath(@__DIR__, "data", "JC_$(JT)_Lx$(Lx)_Ly$(Ly).dat")))
end

# spinglass
function dump_Js(jtype::Val{JT}) where JT
    for i=4:2:32
        J = generate_J(jtype, i*(i-1)*2)
        writedlm(joinpath(@__DIR__, "data", "J_$(JT)_L$i.dat"), J)
    end
end

function load_J(L::Int, jtype::Val{JT}) where JT
    vec(readdlm(joinpath(@__DIR__, "data", "J_$(JT)_L$L.dat")))
end

if abspath(PROGRAM_FILE) == @__FILE__
    dump_JCs(Val(:randn))
    dump_Js(Val(:randn))
end

