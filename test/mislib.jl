using EliminateGraphs
using TropicalTensors
using TropicalTensors: inferier_table
using OMEinsum
using Test

@testset "tropical contract" begin
    t1 = Tropical.(randn(20, 20, 4))
    t2 = Tropical.(randn(4, 20, 30))
    t3 = ein"abc,cbd->abd"(t1,t2)
    @test t3 isa Array{Tropical{Float64},3}
end

@testset "mislib" begin
    T = TMatrix
    g1 = EliminateGraph(5, [1=>2, 2=>3, 2=>4, 3=>4, 4=>5])
    @test mis2(g1) == 3

    # show the final result
    res = ein"ab,bc,cd,bd,de->"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))[].n
    res = ein"(((ab,bc),cd),bd),de->"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))[].n
    @test res == 3

    # show the full table
    tbl = ein"ab,bc,cd,bd,de->abcde"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))

    @test !isinferier(tbl, (1,1,1,1,1), (1,1,1,1,1))
    @test isinferier(tbl, (1,2,2,2,1), (1,1,1,1,1))
    @test isinferier(tbl, (1,2,2,2,1))
    @test !isinferier(tbl, (1,1,1,1,1))

    inferier_table(tbl, (1,1,1,1,1))
    inferier_table(tbl)

    @show tbl
    @test sum(tbl .> Tropical(-555555)) == 12
    @test sum(inferier_table(tbl)) == 32 - 12

    # show the part table, (a, b) are outer legs
    tbl = ein"ab,bc,cd,bd,de->ab"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))

    @test sum(tbl .> Tropical(-555555)) == 3
    @show tbl
    @test isinferier(tbl, (1,2), (1,1))
    @test sum(inferier_table(tbl)) == 2
end

#=
# sub-graph: triangles
# 1.2.1
ein"ij,jk,ki->jk"(T(1//2,1//2), T(1//2,1//2), T(1//2,1//2))

tsplib, satlib, qbflib, junction tree
# 1.2.2.1
ein"ij,ik,jl,kl->l"(T(1//2,1//2), T(1//2,1//2), T(1//2,1//2), T(1//2,1//2))
# 1.2.2.2
ein"ij,ik->jk"(T(1//2,1//1), T(1//2,1//1))
ein"ij,ik,jl,kl,km->lm"(T(1//2,1//2), T(1//2,1//3), T(1//2,1//2), T(1//3,1//2), T(1//3,1//1))
ein"ij,ik,jl,kn,km->lm"(T(1//2,1//2), T(1//2,1//3), T(1//2,1//1), T(1//3,1//1), T(1//3,1//1))
ein"ij,ik,jl,km->lm"(T(1//2,1//2), T(1//2,1//2), T(1//2,1//1), T(1//2,1//1))

# graph K33 + 1, verify the mirror rule: include v, or eliminate v and its mirrors.
T33 = T(1//3, 1//3)
T34 = T(1//3, 1//4)
ein"ia,ib,ic,ja,jb,jc,ka,kb,kc,ab->ijk"(T34, T34, T33,T34,T34,T33,T34,T34,T33,T(1//4,1//4))
=#
