abstract type SiteCollection end
abstract type AbstractLattice <: SiteCollection end
abstract type AbstractSquareLattice <:AbstractLattice end

function Base.getindex(lt::AbstractSquareLattice, i::Int)
    lt[CartesianIndices(size(lt))[i].I...]
end

Base.lastindex(lt::AbstractSquareLattice, i::Int) = size(lt, i)
bravais_size(sq::AbstractSquareLattice) = (sq.Nx, sq.Ny)
function bravais_size(sq::AbstractSquareLattice, i::Int)
    if i==1
        size(sq)[1]
    elseif i==2
        size(sq)[2]
    else
        throw(DimensionMismatch("expected dimension 1 or 2, got $i."))
    end
end
Base.size(sq::AbstractSquareLattice, args...) = bravais_size(sq, args...)

function bond(lt::AbstractSquareLattice, loc1, loc2)
    lt[loc1...], lt[loc2...]
end
