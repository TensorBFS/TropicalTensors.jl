using TropicalYao
using NiLogLikeNumbers
using NiLang, NiLang.AD
using Yao
using Compose

export opt_config, isolve, isolve_largemem

include("square.jl")
include("chimera.jl")

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
