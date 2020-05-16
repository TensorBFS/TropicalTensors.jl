export sgbonds

function square_solve(reg::ArrayReg{B,Tropical{T}}, Lx::Int, Ly::Int, Js::AbstractVector, hs::AbstractVector) where {B,T}
    println("Layer 1/$Lx")
    Js = copy(Js)
    hs = copy(hs)
    for j=1:Ly
        reg |> put(Ly, j=>Gh(T, hs |> popfirst!))
    end
    for j=1:Ly-1
        reg |> put(Ly, (j,j+1)=>G4(T, Js |> popfirst!))
    end
    for i=2:Lx
        println("Layer $i/$Lx")
        for j=1:Ly
            reg |> put(Ly, j=>G2(T, Js |> popfirst!))
        end
        for j=1:Ly
            reg |> put(Ly, j=>Gh(T, hs |> popfirst!))
        end
        for j=1:Ly-1
            reg |> put(Ly, (j,j+1)=>G4(T, Js |> popfirst!))
        end
    end
    sum(state(reg))
end

function solve(sg::Spinglass{LT,T}; usecuda=false) where {LT<:SquareLattice,T}
    # Yao gates
    lt = sg.lattice
    reg = _init_reg(T, lt.Ny, Val(usecuda))
    square_solve(reg, lt.Nx, lt.Ny, sg.Js, sg.hs)
end

function solve(lattice::Viznet.AbstractLattice, Js::Vector{T}, hs::Vector{T}; usecuda=false) where {T}
    sg = Spinglass(lattice, Js, hs)
    solve(sg; usecuda=usecuda)
end

# to index `h`.
function sgvertexorder(lt::SquareLattice)
    v = Int[]
    LI = LinearIndices(lt)
    for i=1:lt.Nx
        for j=1:lt.Ny
            push!(v, LI[i,j])
        end
    end
    return v
end

# to index `J`.
function sgbonds(lt::SquareLattice)
    edges  =Tuple{Int,Int}[]
    append!(edges, v_bonds(lt, 1))
    for i=2:lt.Nx
        append!(edges, h_bonds(lt, i-1))
        append!(edges, v_bonds(lt, i))
    end
    edges
end

function v_bonds(lt::SquareLattice, i::Int)
    LI = LinearIndices(size(lt))
    map(j -> (LI[i,j], LI[i,j+1]), 1:lt.Ny-1)
end

function h_bonds(lt::SquareLattice, i::Int)
    LI = LinearIndices(size(lt))
    map(j -> (LI[i,j], LI[i+1,j]), 1:lt.Ny)
end

regsize(lt::SquareLattice) = lt.Ny
