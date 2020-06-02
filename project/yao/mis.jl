using TropicalTensors
using LightGraphs
using Viznet, Compose
using Viznet: AbstractLattice

struct MISProblem{LT,T} <: AbstractSpinglass{LT,T}
    lattice::LT
    graph::SimpleGraph{Int}
    vertices::Vector{Int}
    bonds::Vector{Tuple{Int,Int}}
    hs::Vector{T}
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

function MISProblem(::Type{T}, lt::AbstractLattice) where T
    MISProblem(lt, zeros(T, length(lt)))
end

function MISProblem(lt::AbstractLattice, hs::AbstractVector{T}) where T
    MISProblem(lt, lattice2graph(lt), sgvertices(lt), sgbonds(lt), hs)
end

function TropicalTensors.bondtensor(::Type{TT}, mis::MISProblem, i::Int) where TT
    a, b = mis.bonds[i]
    ia = findfirst(==(a), mis.vertices)
    ib = findfirst(==(b), mis.vertices)
    [one(TT) TT(1/degree(mis.graph, ia)); TT(1/degree(mis.graph, ib)) zero(TT)]
end

function TropicalTensors.vertextensor(::Type{TT}, mis::MISProblem, i::Int) where TT
    h = mis.hs[i]
    TT.([h, -h])
end

using Test, Random
@testset "lattice to misproblem" begin
    lt = SquareLattice(5, 5)
    g = lattice2graph(lt)
    @test ne(g) == 40
    mis = MISProblem(Float64, lt)
    @test mis.graph == g
    @test bondtensor(Tropical{Float64}, mis, 3) isa Matrix{Tropical{Float64}}
end

@testset "spinglass" begin
    lt = SquareLattice(10, 8)
    sg = MISProblem(Float64, lt)
    res = solve(sg; usecuda=false)
    @test res.n ≈ 40
end

@testset "spinglass" begin
    Random.seed!(4)
    lt = rand_maskedsquare(10, 8, 1.0)
    sg = MISProblem(Float64, lt)
    res = solve(sg; usecuda=false)
    @test res.n ≈ 20
end

struct MISOptConfig{LT,T}
    mis::MISProblem{LT}
    eng::T
    grad_hs::Vector{T}
end

function vizgrad_mis(mis::MISProblem, grad_hs::AbstractVector; r=0.015)
    lt = mis.lattice
    nb1 = compose(nodestyle(:default; r=r), fill("white"), stroke("black"), linewidth(0.4mm))
    nb2 = compose(nodestyle(:default; r=r), fill("black"), stroke("black"), linewidth(0.4mm))
    eb1 = compose(bondstyle(:default), linewidth(0.7mm), stroke("skyblue"))
    cdots = canvas() do
        for i in sgvertices(lt)
            ii = findfirst(==(i), mis.vertices)
            if grad_hs[ii] > 0
                nb1 >> lt[i]
            elseif grad_hs[ii] < 0
                nb2 >> lt[i]
            else
                error("index $i not set!")
            end
        end
        for (i,j) in sgbonds(lt)
            eb1 >> lt[i;j]
        end
    end
    compose(context(), cdots)
end

function Base.display(sgres::MISOptConfig)
    Base.display(vizgrad_mis(sgres.mis, sgres.grad_hs))
end

using ForwardDiff
function opt_config_forwarddiff(mis::MISProblem{LT, T}; usecuda=false) where {LT, T}
    ghs = ForwardDiff.gradient(hs->solve(MISProblem(mis.lattice, hs), usecuda=usecuda).n, mis.hs)
    MISOptConfig(mis, solve(mis, usecuda=usecuda).n, ghs)
end

Random.seed!(4)
lt = rand_maskedsquare(15, 15, 0.5)
mis = MISProblem(Float64, lt)
res = opt_config_forwarddiff(mis)
