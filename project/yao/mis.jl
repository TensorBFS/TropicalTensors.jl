using TropicalTensors
using LightGraphs
using Viznet
using Viznet: AbstractLattice

struct MISProblem{LT} <: AbstractSpinglass{LT}
    lattice::LT
    graph::SimpleGraph{Int}
    vertices::Vector{Int}
    bonds::Vector{Tuple{Int,Int}}
end

function lattice2graph(lt::AbstractLattice)
    g = SimpleGraph(length(lt))
    vertices = sgvertices(lt)
    for (i, j) in sgbonds(lt)
        ia = findfirst(==(i), vertices)
        ib = findfirst(==(j), vertices)
        add_edge!(g, ia, ib)
    end
    g
end

MISProblem(lt::AbstractLattice) = MISProblem(lt, lattice2graph(lt), sgvertices(lt), sgbonds(lt))

function TropicalTensors.bondtensor(::Type{TT}, mis::MISProblem, i::Int) where TT
    a, b = mis.bonds[i]
    ia = findfirst(==(a), mis.vertices)
    ib = findfirst(==(b), mis.vertices)
    [one(TT) TT(1/degree(mis.graph, ia)); TT(1/degree(mis.graph, ib)) zero(TT)]
end
function TropicalTensors.vertextensor(::Type{TT}, mis::MISProblem, i::Int) where TT
    ones(TT, 2)
end

function TropicalTensors.solve(sg::MISProblem{LT}; usecuda=false) where {LT}
    solve(Tropical{Float64}, sg; usecuda=usecuda)
end

using Test, Random
@testset "lattice to misproblem" begin
    lt = SquareLattice(5, 5)
    g = lattice2graph(lt)
    @test ne(g) == 40
    mis = MISProblem(lt)
    @test mis.graph == g
    @test bondtensor(Tropical{Float64}, mis, 3) isa Matrix{Tropical{Float64}}
end

@testset "spinglass" begin
    lt = SquareLattice(10, 8)
    sg = MISProblem(lt)
    res = solve(sg; usecuda=false)
    @test res.n ≈ 40
end

@testset "spinglass" begin
    Random.seed!(4)
    lt = rand_maskedsquare(10, 8, 1.0)
    sg = MISProblem(lt)
    res = solve(sg; usecuda=false)
    @show res.n
    @test res.n ≈ 20
end
