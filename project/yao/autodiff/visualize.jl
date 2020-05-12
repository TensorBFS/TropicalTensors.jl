using Compose
using Viznet

set_default_graphic_size(20cm, 20cm)

function spinglass(lt, coulping::Vector, grid)
    nb1 = compose(nodestyle(:default; r=0.008), fill("white"), stroke("black"), linewidth(0.4mm))
    nb2 = compose(nodestyle(:default; r=0.008), fill("black"), stroke("white"), linewidth(0.4mm))
    eb1 = compose(bondstyle(:default), linewidth(0.7mm), stroke("skyblue"))
    eb2 = compose(bondstyle(:default), linewidth(0.7mm), stroke("orange"))
    _coupling = copy(coulping)
    eb() = popfirst!(_coupling) > 0 ? eb1 : eb2
    g = canvas() do
        for j=1:lt.Ny
            for i=1:lt.Nx
                (grid[i,j] > 0 ? nb1 : nb2) >> lt[i,j]
            end
        end
        for i=1:lt.Nx-1
            eb() >> lt[(i,1); (i+1,1)]
        end
        for j=2:lt.Ny
            for i=1:lt.Nx
                eb() >> lt[(i,j-1); (i,j)]
            end
            for i=1:lt.Nx-1
                eb() >> lt[(i,j); (i+1,j)]
            end
        end
    end
    compose(context(), g)
end

function assign_grid(L, g::AbstractVector)
    grid = zeros(Int, L, L)
    grid[1,1] = 1
    println("Layer 1/$L")
    k = 0
    for i=1:L-1
        k += identity(1)
        assign_one!(grid, (i,1), (i+1,1), g[k])
    end
    for j=2:L
        println("Layer $j/$L")
        for i=1:L
            k += identity(1)
            assign_one!(grid, (i,j-1), (i,j), g[k])
        end
        for i=1:L-1
            k += identity(1)
            assign_one!(grid, (i,j), (i+1,j), g[k])
        end
    end
    return grid
end

function assign_one!(grid, x, y, g)
    if grid[x...] == 0 && grid[y...] == 0
        error("")
    elseif grid[x...] == 0
        grid[x...] = sign(g)*grid[y...]
    elseif grid[y...] == 0
        grid[y...] = sign(g)*grid[x...]
    else
        @assert grid[y...] == sign(g)*grid[x...]
    end
end

include("../datadump.jl")
function show_coupling(::Type{T}, L; jtype, filename="_lattice.svg") where T
    Js = load_J(L, Val(jtype))
    lt = SquareLattice(L, L)
    g = readdlm(joinpath(@__DIR__, "../data", "spinglassgrad_$(jtype)_L$(L)_$(T).dat"))
    grid = Int.((assign_grid(L, vec(g)) .+ 1) ./ 2)
   spinglass(lt, Js, grid)
end

show_coupling(Int32, 28; jtype=:randpm)
