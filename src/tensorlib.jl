export hypercubicI, spinglass_bond_tensor, spinglass_vertex_tensor, spinglass_g4_tensor

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

function spinglass_bond_tensor(Jij::Real)
    tij = Tropical(Jij)
    _tij = Tropical(-Jij)
    return [tij _tij; _tij tij]
end

function spinglass_g4_tensor(Jij::Real)
    reshape(ein"ab->abab"(spinglass_bond_tensor(Jij)), 4, 4)
end
