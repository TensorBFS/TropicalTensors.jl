export spinglass_g4_tensor!, spinglass_bond_tensor!, spinglass_g16_tensor!

"""
    spinglass_bond_tensor!(mat, Jij)

`mat` should be a `one` tensor.
"""
@i function spinglass_bond_tensor!(mat::AbstractMatrix, Jij::Real)
    @safe @assert size(mat) == (2,2)
    Tropical(Jij)
    muleq(mat[1,1], Jij)
    muleq(mat[2,2], Jij)
    (~muleq)(mat[1,2], Jij)
    (~muleq)(mat[2,1], Jij)
    (~Tropical)(Jij)
end

@i function spinglass_g4_tensor!(mat::Diagonal, Jij::Real)
    @safe @assert size(mat) == (4,4)
    Tropical(Jij)
    muleq(mat[1,1], Jij)
    (~muleq)(mat[2,2], Jij)
    (~muleq)(mat[3,3], Jij)
    muleq(mat[4,4], Jij)
    (~Tropical)(Jij)
end

@i function spinglass_g16_tensor!(out!::AbstractMatrix{T}, Js) where T<:Tropical
    @safe @assert length(Js) == 16
    @safe @assert size(out!) == (16, 16)
    @routine begin
        y ← ones(T,2,2,2,2,2,2,2,2)
        xs ← ([ones(T,2,2) for i=1:16]...,)
        for i = 1:length(Js)
            spinglass_bond_tensor!(tget(xs,i), Js[i])
        end
        ieinsum!(ein"aα,aβ,aγ,aδ,bα,bβ,bγ,bδ,cα,cβ,cγ,cδ,dα,dβ,dγ,dδ->abcdαβγδ", xs, y)
    end
    for i=1:length(out!)
        muleq(out![i], y[i])
    end
    ~@routine
end
