include("spinglass.jl")

L = 28

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

@time solve_and_dump(Int32, L; jtype=:randpm)
#analyze_grad(Int32, L; jtype=:randpm)

using Viznet

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

show_coupling(Int16, L; jtype=:randpm)
