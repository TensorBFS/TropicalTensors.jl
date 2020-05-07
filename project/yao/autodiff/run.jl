include("spinglass.jl")

function solve_and_dump(::Type{T}, L; jtype) where T
    eng, g = opt_config(T, L; jtype=jtype)
    filename = joinpath(@__DIR__, "../data", "spinglassgrad_$(jtype)_L$(L)_$(T).dat")
    println("optimal energy is $eng, saving gradients to $filename")
    writedlm(filename, g)
end

#function analyze_grad(::Type{T}, L; jtype) where T
#    g = readdlm(joinpath(@__DIR__, "../data", "spinglassgrad_$(jtype)_L$(L)_$(T).dat"))
#    grid = assign_grid(L, vec(g))
#    print_grid(grid)
#end

solve_and_dump(Int32, 10; jtype=:randpm)
#analyze_grad(Int32, 10; jtype=:randpm)

include("../../visualizetools/compose.jl")

function show_coupling(::Type{T}, L; jtype, filename="_lattice.svg") where T
    Js = load_J(L, Val(jtype))
    lt = SquareLattice(L, L, 0.2)
    bonds = get_bonds(lt.Nx, lt.Ny)
    g = readdlm(joinpath(@__DIR__, "../data", "spinglassgrad_$(jtype)_L$(L)_$(T).dat"))
    grid = Int.((assign_grid(L, vec(g)) .+ 1) ./ 2)
    pbonds = bonds[Js .> 0]
    nbonds = bonds[Js .< 0]
    showbonds(10, 10, grid, [("red"=>pbonds), "green"=>nbonds])
end

function get_bonds(Lx, Ly)
    bonds = Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]
    for i=1:Lx-1
        push!(bonds, ((i,1), (i+1,1)))
    end
    for j=2:Ly
        for i=1:Lx
            push!(bonds, ((i,j-1), (i,j)))
        end
        for i=1:Lx-1
            push!(bonds, ((i,j), (i+1,j)))
        end
    end
    return bonds
end

show_coupling(Int32, 10; jtype=:randpm)
