using Compose

set_default_graphic_size(4cm, 4cm)

struct SquareLattice
    Nx::Int
    Ny::Int
    r::Float64
end

function Base.getindex(lt::SquareLattice, i::Int, j::Int)
    step = unit(lt)
    (i-0.5 + lt.r)*step, (j-0.5 + lt.r)*step
end

function Base.getindex(lt::SquareLattice, i::Int)
    lt[CartesianIndices((lt.Nx,lt.Ny))[i].I...]
end
Base.lastindex(lt::SquareLattice, i::Int) = size(lt,i)
Base.size(sq::SquareLattice) = (sq.Nx, sq.Ny)
Base.size(sq::SquareLattice, i::Int) = i==1 ? sq.Nx : sq.Ny
function bond(lt::SquareLattice, loc1, loc2)
    lt[loc1...], lt[loc2...]
end

unit(lt::SquareLattice) = 1/(max(lt.Nx, lt.Ny)+2*lt.r)

using Test
@testset "sq lattice" begin
    lt = SquareLattice(10, 20, 0.2)
    @test size(lt, 1) == 10
    @test size(lt, 2) == 20
    @test size(lt) == (10, 20)
    @test all(lt[1,1] .≈ (unit(lt) * (lt.r+0.5), unit(lt) * (lt.r+0.5)))
    lt = SquareLattice(10, 10, 0.2)
    @test all(lt[end,end] .≈ (1-unit(lt) * (lt.r+0.5), 1-unit(lt) * (lt.r+0.5)))
end

function locs(sq::SquareLattice)
    xs = Float64[]
    ys = Float64[]
    for j=1:size(sq, 2), i=1:size(sq, 1)
        xi, yi = sq[i,j]
        push!(xs, xi)
        push!(ys, yi)
    end
    return xs, ys
end

rs(lt::SquareLattice) = fill(lt.r*unit(lt), lt.Nx, lt.Ny)

function lattice(lt)
    circle(locs(lt)..., rs(lt))
end

using Colors
function showbonds(Nx::Int, Ny::Int, color_sites::AbstractMatrix, color_bonds::AbstractVector;
                filename="_lattice.svg")
    lt = SquareLattice(Nx, Ny, 0.3)
    scolors = LCHab.(vec(color_sites).*200, 230, 57)
    composition = compose(context(),
        (context(), lattice(lt), fill(scolors), stroke("silver"), linewidth(0.05mm)),
        #(context(), circle(lt[1,1]..., 0.04)),
        #(context(), circle(lt[2,1]..., 0.04)),
        map(cb->(context(), line(map(x->bond(lt, x...), cb.second)), stroke(cb.first), linewidth(0.5mm)), color_bonds)...
        )
    composition |> SVG(filename)
end

# showbonds(10, 10, [("red"=>[(2,3), (4,5)])])
