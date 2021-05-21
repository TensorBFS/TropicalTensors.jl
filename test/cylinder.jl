using TropicalTensors
using ForwardDiff
using Test

@testset "spinglass" begin
    lt = Cylinder(10, 8)
    sg = Spinglass(lt, ones(Float32, 152), zeros(Float32, 80))
    res = solve(sg; usecuda=false)
    @test res.n == 152
end

@testset "counting tropical" begin
    lt = Cylinder(10, 8)
    sg = Spinglass(lt, ones(Float32, 152), zeros(Float32, 80))
    res = solve(CountingTropicalF32, sg; usecuda=false)
    @test res.n == 152
    @test res.c == 2
end

using BitBasis
function count_degeneracy_exact(sg::Spinglass{LT,T}) where {LT, T}
    elist = T[]
    bonds = sgbonds(sg.lattice)
    for b in basis(BitStr64{length(sg.lattice)})
        eng = zero(T)
        for ((i, j), J) in zip(bonds, sg.Js)
            eng += b[i] == b[j] ? J : -J
        end
        push!(elist, eng)
    end
    findall(==(minimum(elist)), elist) |> length
end

@testset "counting tropical 2" begin
    lt = Cylinder(4, 4)
    sg = Spinglass(lt, rand([-1, 1], 28), zeros(Int, 16))
    res = solve(CountingTropical{Int,Int}, sg; usecuda=false)
    @test res.c == count_degeneracy_exact(sg)
end

@testset "forwarddiff" begin
    L = 3
    Js = ones(Float32, 2*L^2-L)
    gs = ForwardDiff.gradient(x->solve(Cylinder(L, L), x, zeros(eltype(x), L^2); usecuda=false).n, Js)
    @test gs ≈ ones(Float32, 2*L^2-L)
end
