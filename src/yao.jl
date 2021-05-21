using Yao
using LuxurySparse
using LinearAlgebra
using TropicalYao
using Viznet, Compose

export solve, SquareLattice, ChimeraLattice
export sgbonds, sgvertices

export Gh, Gvb, Ghb, G16, Gcut, Gcp, Greset
export vertextensor, bondtensor

function _init_reg(::Type{T}, L::Int, usecuda::Val{:false}) where T
    ArrayReg(ones(T, 1<<L))
end

function solve(sg::AbstractSpinglass{LT,T}; usecuda=false) where {LT,T}
    solve(Tropical{T}, sg; usecuda=usecuda)
end

function solve_and_count(sg::Spinglass{LT,T}; usecuda=false) where {LT,T}
    solve(CountingTropical{T,T}, sg; usecuda=usecuda)
end

struct SpinglassOptConfig{LT,T}
    sg::Spinglass{LT,T}
    eng::T
    grad_Js::Vector{T}
    grad_hs::Vector{T}
end

export tofile, fromfile
function fromfile(::Type{T}, filename::String, lattice) where T
    nv = length(lattice)
    ne = length(sgbonds(lattice))
    Js = zeros(T, ne)
    hs = zeros(T, nv)
    Js_grad = zeros(T, ne)
    hs_grad = zeros(T, nv)
    open(filename, "r") do f
        for i=1:length(Js)
            l = strip(readline(f))
            j, g = split(l, ' ')
            Js[i] = parse(T, j)
            Js_grad[i] = parse(T, g)
        end
        for i=1:length(hs)
            l = strip(readline(f))
            h, g = split(l, ' ')
            hs[i] = parse(T, h)
            hs_grad[i] = parse(T, g)
        end
    end
    SpinglassOptConfig(Spinglass(lattice, Js, hs), T(0), Js_grad, hs_grad)
end

function tofile(filename::String, res::SpinglassOptConfig)
    sg = res.sg
    open(filename, "w") do f
        for (J, g) in zip(sg.Js, res.grad_Js)
            write(f, "$J $g\n")
        end
        for (h, g) in zip(sg.hs, res.grad_hs)
            write(f, "$h $g\n")
        end
    end
    return res
end

tofile(filename::String) = res->tofile(filename, res)

include("square.jl")
include("chimera.jl")
include("second_neighbor.jl")
include("cylinder.jl")
include("cubic.jl")

function assign_Js(lt::Viznet.AbstractLattice, g::AbstractVector{T}) where T
    grid = zeros(length(lt))
    grid[1] = 1
    configs = length(lt)
    bonds = sgbonds(lt)
    remain_bonds = collect(1:length(bonds))
    for i=1:10
        nrem = length(remain_bonds)
        for (i_, k) in enumerate(Base.Iterators.reverse(copy(remain_bonds)))
            i, j = bonds[k]
            res = assign_one!(grid, i, j, g[k])
            if res
                deleteat!(remain_bonds, nrem-i_+1)
            end
        end
        isempty(remain_bonds) && break
    end
    if !isempty(remain_bonds)
        error("bonds $remain_bonds can not be set!")
    end
    return grid
end

function assign_Js_hs(lt::Viznet.AbstractLattice, grad_Js::AbstractVector{T}, grad_hs) where T
    grid = zeros(length(lt))
    for (v, hg) in zip(sgvertices(lt), grad_hs)
        grid[v] = hg
    end
    for ((i, j), Jg) in zip(sgbonds(lt), grad_Js)
        assign_one!(grid, i, j, Jg)
    end
    return grid
end

function assign_one!(grid, x, y, g)
    if grid[x] == 0 && grid[y] == 0
        return false
    elseif grid[x] == 0
        grid[x] = sign(g)*grid[y]
    elseif grid[y] == 0
        grid[y] = sign(g)*grid[x]
    else
        if grid[y] != sign(g)*grid[x]
            @warn "Grid assign σ($y) = $(grid[y]) and σ($x) = $(grid[x]) inconsistent with bond gradient $(g)!"
        end
    end
    return true
end
