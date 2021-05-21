using TropicalTensors
using TropicalTensors: getbonds, getvertices, rand_Js, rand_hs
using Test
using Random

@testset "spinglass" begin
    lt = CubicLattice(4,4,4)
    sg = Spinglass(lt, rand_Js(Float32, Ferro(), lt), rand_hs(Float32, Zero(), lt))
    res = solve(sg; usecuda=false)
    @test res.n == 144
end

@testset "counting tropical" begin
    lt = CubicLattice(4,4,4)
    sg = Spinglass(lt, rand_Js(Float32, Ferro(), lt), rand_hs(Float32, Zero(), lt))
    res = solve(CountingTropicalF32, sg; usecuda=false)
    @test res.n == 144
    @test res.c == 2
end

using BitBasis
function count_degeneracy_exact(sg::Spinglass{LT,T}) where {LT, T}
    elist = T[]
    bonds = getbonds(sg.lattice)
    lt = sg.lattice
    LI = LinearIndices((lt.Nx, lt.Ny, lt.Nz))
    for b in basis(BitStr64{length(sg.lattice)})
        eng = zero(T)
        for ((i, j), J) in sg.Js
            eng += (b[LI[i...]] == b[LI[j...]]) ? J : -J
        end
        push!(elist, eng)
    end
    findall(==(minimum(elist)), elist) |> length
end

@testset "counting tropical 2" begin
    Random.seed!(3)
    lt = CubicLattice(2, 3, 3)
    sg = Spinglass(lt, rand_Js(Float32, Randpm(), lt), rand_hs(Float32, Zero(), lt))
    res = solve(CountingTropical{Int,Int}, sg; usecuda=false)
    @test res.c == count_degeneracy_exact(sg)
end