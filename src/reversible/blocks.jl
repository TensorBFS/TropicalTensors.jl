export spinglass_g4_tensor!, spinglass_bond_tensor!, spinglass_g16_tensor!

"""
    spinglass_bond_tensor!(mat, Jij)

`mat` should be a `one` tensor.
"""
@i function spinglass_bond_tensor!(mat::AbstractMatrix, Jij::Real)
    @safe @assert size(mat) == (2,2)
    Tropical(Jij)
    mat[1,1] *= identity(Jij)
    mat[2,2] *= identity(Jij)
    mat[1,2] /= identity(Jij)
    mat[2,1] /= identity(Jij)
    (~Tropical)(Jij)
end

@i function spinglass_g4_tensor!(mat::Diagonal, Jij::Real)
    @safe @assert size(mat) == (4,4)
    Tropical(Jij)
    mat[1,1] *= identity(Jij)
    mat[2,2] /= identity(Jij)
    mat[3,3] /= identity(Jij)
    mat[4,4] *= identity(Jij)
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
        ixs ← ([(ix...,) for ix in split("aα,aβ,aγ,aδ,bα,bβ,bγ,bδ,cα,cβ,cγ,cδ,dα,dβ,dγ,dδ", ',')]...,)
        iy ← ("abcdαβγδ"...,)
        NiLog.einsum!(ixs, xs, iy, y)
    end
    for i=1:length(out!)
        out![i] *= identity(y[i])
    end
    ~@routine
end
