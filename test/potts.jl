using TropicalTensors
using TropicalTensors: potts_randpm_J, SquareLattice
using Test, Random

function exact_solve_potts3(lt::SquareLattice, jd::Dict)
    eng_list = Float64[]
    for ci in CartesianIndices((fill(3, lt.Nx*lt.Ny)...,))
        push!(eng_list, potts3_energy(lt, jd, ci.I))
    end
    max_energy = maximum(eng_list)
    degeneracy = count(==(max_energy), eng_list)
    return max_energy, degeneracy
end

function potts3_energy(lt::SquareLattice, jd::Dict, cfg)
    LI = LinearIndices((lt.Nx, lt.Ny))
    eng = 0.0
    for i=1:lt.Nx
        for j=1:lt.Ny
            i!=lt.Nx && if cfg[LI[i,j]] == cfg[LI[i+1,j]]
                eng += jd[(i,j)=>(i+1,j)]
            else
                eng -= 0.5*jd[(i,j)=>(i+1,j)]
            end
            j!=lt.Ny && if cfg[LI[i,j]] == cfg[LI[i,j+1]]
                eng += jd[(i,j)=>(i,j+1)]
            else
                eng -= 0.5*jd[(i,j)=>(i,j+1)]
            end
        end
    end
    return eng
end

@testset "check potts" begin
    L = 9
    lt = SquareLattice(L, L)
    CI = CartesianIndices(lt)
    for i=1:10
        Random.seed!(i)
        dj = potts_randpm_J(lt)
        sj = Float64[]
        for b in sgbonds(lt)
            push!(sj, dj[CI[b[1]].I => CI[b[2]].I])
        end
        res1 = solve_potts(CountingTropical{Float64}, Val(2), lt, dj; usecuda=false)
        res2 = solve(CountingTropical{Float64}, Spinglass(lt, sj, zeros(L*L)); usecuda=false)
        @test res1.n ≈ res2.n
        @test res1.c ≈ res2.c
    end
end

@testset "compare exact" begin
    L = 3
    lt = SquareLattice(L, L)
    jd = potts_randpm_J(lt)
    res1 = solve_potts(CountingTropical{Float64}, Val(3), lt, jd; usecuda=false)
    eng, deg = exact_solve_potts3(lt, jd)
    @test res1.n == eng
    @test res1.c == deg
end