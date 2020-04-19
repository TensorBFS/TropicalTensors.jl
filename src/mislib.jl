export TMatrix, isinferier

TMatrix(a::T, b::T) where T = [Tropical(T(0)) Tropical(b); Tropical(a) Tropical(T(TROPICAL_ZERO))]

isinferier(tbl, σA::NTuple{X,Int}, σB::NTuple{X,Int}) where X = σA!=σB && tbl[σA...].n <= tbl[σB...].n && all(σA .>= σB)
isinferier(tbl, σA::NTuple{X,Int}) where X = isempty(σA) || any(σB -> isinferier(tbl, σA, σB.I), CartesianIndices(tbl))

# this marks eliminatable elements
inferier_table(tbl, σA::NTuple{X,Int}) where X = map(σB -> isinferier(tbl, σA, σB.I), CartesianIndices(tbl))
inferier_table(tbl) where X = map(σB -> isinferier(tbl, σB.I), CartesianIndices(tbl))
