abstract type Jtype end

export Spinglass, Randn, Rand, Randpm, Ferro, Zero, rand_spinglass

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

struct Spinglass{LT,JT}
    lattice::LT
    Js::Vector{JT}
    hs::Vector{JT}
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
