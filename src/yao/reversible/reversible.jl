using TropicalYao
using NiLogLikeNumbers
using NiLang, NiLang.AD
using Yao
using Compose

export opt_config, isolve, isolve_largemem

include("square.jl")
include("chimera.jl")

struct SpinglassOptConfig{LT,T}
    sg::Spinglass{LT,T}
    eng::T
    grad_J::Vector{T}
    grad_h::Vector{T}
end

regsize(lt::SquareLattice) = lt.Ny
regsize(lt::ChimeraLattice) = lt.Ny*4

function opt_config(sg::Spinglass{LT,T}) where {LT,T}
    nbit = regsize(sg.lattice)
    reg = ArrayReg(ones(Tropical{T}, 1<<nbit))
    A = stack4reg(reg, nbit)
    B = stack4reg(reg, sg.lattice.Nx-1)
    eng, sg, reg, A, B = isolve(T(0.0), sg, reg, A, B)
    sgg = Spinglass(sg.lattice, GVar.(sg.Js, zero(sg.Js)), GVar.(sg.hs, zero(sg.hs)))
    gres = (~isolve)(GVar(eng, T(1)), sgg, GVar(reg), GVar(A), GVar(B))
    empty!(NiLang.GLOBAL_STACK)
    return SpinglassOptConfig(sg, eng, grad.(gres[2].Js), grad.(gres[2].hs))
end

function Base.display(sgres::SpinglassOptConfig)
    Base.display(vizgrad_J(sgres.sg, sgres.grad_J))
end

export vizgrad_J
function vizgrad_J(sg::Spinglass, grad_J::AbstractVector)
    lt = sg.lattice
    grid = assign_grid(lt, grad_J)
    nb1 = compose(nodestyle(:default; r=0.015), fill("white"), stroke("black"), linewidth(0.4mm))
    nb2 = compose(nodestyle(:default; r=0.015), fill("black"), stroke("white"), linewidth(0.4mm))
    eb1 = compose(bondstyle(:default), linewidth(0.7mm), stroke("skyblue"))
    eb2 = compose(bondstyle(:default), linewidth(0.7mm), stroke("orange"))
    cdots = canvas() do
        for i=1:length(lt)
            if grid[i] > 0
                nb1 >> lt[i]
            elseif grid[i] < 0
                nb2 >> lt[i]
            end
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
