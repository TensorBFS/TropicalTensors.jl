abstract type Jtype end

export AbstractSpinglass, Spinglass, Randn, Rand, Randpm, Ferro, Zero, rand_spinglass

struct Randn <: Jtype end
struct Rand <: Jtype end
struct Randpm <: Jtype end
struct Ferro <: Jtype end
struct Zero <: Jtype end

_get_J(::Ferro) = 1.0
_get_J(::Randn) = randn()
_get_J(::Rand) = rand()
_get_J(::Randpm) = rand([-1,1])
_get_J(::Zero) = 0.0

abstract type AbstractSpinglass{LT, T} end

struct Spinglass{LT,T} <: AbstractSpinglass{LT,T}
    lattice::LT
    Js::Vector{T}
    hs::Vector{T}
end

function bondtensor(::Type{TT}, sg::Spinglass, i::Int) where TT
    J = sg.Js[i]
    TT.([J -J; -J J])
end
function vertextensor(::Type{TT}, sg::Spinglass, i::Int) where TT
    h = sg.hs[i]
    TT.([h, -h])
end

function rand_spinglass(::Type{T}, lt::Viznet.AbstractLattice;
                    jt::Jtype=Randn(), ht::Jtype=Zero(), seed::Int=2,
                    dumpto=nothing) where T
    Random.seed!(seed)
    hs = T[_get_J(ht) for i=1:length(lt)]
    Js = T[_get_J(jt) for i=1:length(sgbonds(lt))]
    sg = Spinglass(lt, Js, hs)
    if !(dumpto isa Nothing)
        dump_params(dumpto, sg, jt, ht, seed)
    end
    return sg
end

export viz_sg
function viz_sg(sg::Spinglass; r=0.015)
    lt = sg.lattice
    nb1 = compose(nodestyle(:default; r=r), fill("white"), stroke("black"), linewidth(0.4mm))
    eb1 = compose(bondstyle(:default), linewidth(0.7mm), stroke("skyblue"))
    eb2 = compose(bondstyle(:default), linewidth(0.7mm), stroke("orange"))
    tb = textstyle(:default)
    cdots = canvas() do
        for v in sgvertices(lt)
            nb1 >> lt[v]
            tb >> (lt[v], "$v")
        end
        for ((i,j),v) in zip(sgbonds(lt), sg.Js)
            if v > 0
                eb1 >> lt[i;j]
            else
                eb2 >> lt[i;j]
            end
        end
    end
    compose(context(), cdots)
end

function Base.display(sg::Spinglass; r=0.02)
    Base.display(viz_sg(sg; r=r))
end

export eval_onconfig
function eval_onconfig(sg::Spinglass{LT,T}, grid) where {LT,T}
    out = T(0)
    for ((i,j), J) in zip(sg.lattice |> sgbonds, sg.Js)
        out += J * grid[i] * grid[j]
    end
    return out
end
