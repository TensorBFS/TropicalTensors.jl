using Yao
using LuxurySparse
using LinearAlgebra
using TropicalYao
using Viznet

export solve, SquareLattice, ChimeraLattice
export sgbonds

Gh(::Type{T}, h) where T = matblock(Diagonal(spinglass_mag_tensor(T(h))) |> LuxurySparse.staticize)
G2(::Type{T}, J) where T = matblock(spinglass_bond_tensor(T(J)) |> LuxurySparse.staticize)
G4(::Type{T}, J) where T = matblock(Diagonal(spinglass_g4_tensor(T(J))) |> LuxurySparse.staticize)
G16(::Type{T}, Js) where T = matblock(spinglass_g16_tensor(T.(Js)) |> LuxurySparse.staticize)

function _init_reg(::Type{T}, L::Int, usecuda::Val{:false}) where T
    reg = ArrayReg(ones(Tropical{T}, 1<<L))
end

struct SpinglassOptConfig{LT,T}
    sg::Spinglass{LT,T}
    eng::T
    grad_Js::Vector{T}
    grad_hs::Vector{T}
end

include("square.jl")
include("chimera.jl")

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
    vorder = sgvertexorder(lt)
    for (i, hg) in enumerate(grad_hs)
        grid[vorder[i]] = hg
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
        @assert grid[y] == sign(g)*grid[x]
    end
    return true
end

function Base.display(sgres::SpinglassOptConfig)
    Base.display(vizgrad_J(sgres.sg, sgres.grad_Js, sgres.grad_hs))
end

export vizgrad_J
function vizgrad_J(sg::Spinglass, grad_Js::AbstractVector, grad_hs::AbstractVector)
    lt = sg.lattice
    grid = assign_Js_hs(lt, grad_Js, grad_hs)
    nb1 = compose(nodestyle(:default; r=0.015), fill("white"), stroke("black"), linewidth(0.4mm))
    nb2 = compose(nodestyle(:default; r=0.015), fill("black"), stroke("white"), linewidth(0.4mm))
    eb1 = compose(bondstyle(:default), linewidth(0.7mm), stroke("skyblue"))
    eb2 = compose(bondstyle(:default), linewidth(0.7mm), stroke("orange"))
    cdots = canvas() do
        for i=1:length(lt)
            if grid[i] > 0
                nb1 >> lt[i]
            elseif grid[i] < 0
                nb2 >> lt[i]
            else
                error("index $i not set!")
            end
        end
        for ((i,j),v) in zip(sgbonds(lt), sg.Js)
            if v > 0
                eb1 >> lt[i;j]
            else
                eb2 >> lt[i;j]
            end
        end
    end
    compose(context(), cdots)
end
