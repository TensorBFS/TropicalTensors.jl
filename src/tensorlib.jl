export hypercubicI, spinglass_bond_tensor, spinglass_vertex_tensor, spinglass_g4_tensor
export spinglass_g16_tensor

"""
    hypercubicI([T], ndim::Int, D::Int)

Get a `ndim`-dimensional identity hypercubic, with bound dimension `D`.
"""
function hypercubicI(::Type{T}, ndim::Int, D::Int) where T
    res = zeros(T, fill(D, ndim)...)
    @inbounds for i=1:D
        res[fill(i,ndim)...] = one(T)
    end
    return res
end

hypercubicI(ndim::Int, D::Int) = hypercubicI(Float64, ndim, D)

function spinglass_vertex_tensor(::Type{T}, ndim::Int) where T
    hypercubicI(Tropical{T}, ndim, 2)
end

spinglass_vertex_tensor(ndim::Int) = spinglass_vertex_tensor(Float64, ndim)

function spinglass_bond_tensor(Jij)
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

function spinglass_g4_tensor(Jij)
    reshape(ein"ab->abab"(spinglass_bond_tensor(Jij)), 4, 4)
end

function spinglass_g16_tensor(Js)
    @assert length(Js) == 16
    xs = map(spinglass_bond_tensor, Js)
    reshape(ein"(aα,aβ,aγ,aδ),(bα,bβ,bγ,bδ),(cα,cβ,cγ,cδ),(dα,dβ,dγ,dδ)->abcdαβγδ"(xs...), 16, 16)
end
