using DelimitedFiles

_get_J(::Val{:ferro}) = 1.0
_get_J(::Val{:randn}) = randn()
_get_J(::Val{:rand}) = rand()
_get_J(::Val{:randpm}) = rand([-1,1])

function generate_J(type::Val, size...)
    J = zeros(size...)
    for i=1:length(J)
        J[i] = _get_J(type)
    end
    return J
end

# chimera
function dump_JCs(jtype::Val{JT}) where JT
    for i=3:8
        J = generate_J(jtype, i*(i-1)*8 + i^2*16)
        writedlm(joinpath(@__DIR__, "data", "JC_$(JT)_Lx$(i)_Ly$(i).dat"), J)
    end
end

function load_JC(Lx::Int, Ly::Int, jtype::Val{JT}) where JT
    vec(readdlm(joinpath(@__DIR__, "data", "JC_$(JT)_Lx$(Lx)_Ly$(Ly).dat")))
end

# spinglass
function dump_Js(jtype::Val{JT}) where JT
    for i=4:2:32
        J = generate_J(jtype, i*(i-1)*2)
        writedlm(joinpath(@__DIR__, "data", "J_$(JT)_L$i.dat"), J)
    end
end

function load_J(L::Int, jtype::Val{JT}) where JT
    vec(readdlm(joinpath(@__DIR__, "data", "J_$(JT)_L$L.dat")))
end

function load_JC(Lx::Int, Ly::Int, jtype::Val{:ferro})
    return ones(Lx*(Ly-1)*4 + Ly*(Lx-1)*4 + Lx*Ly*16)
end

index0(ix, iy, Lx) = 8(ix-1) + 8Lx*(iy-1)

function red_bond!(bonds, ix::Int, Lx::Int, Ly::Int, Js)
    k = 0
    for i=1:Ly-1
        i0 = index0(ix, i, Lx)
        j0 = index0(ix, i+1, Lx)
        for j = 1:4
            push!(bonds, (i0+j, j0+j))
            k += 1
        end
    end
    for i=1:Ly
        i0 = index0(ix, i, Lx)
        push16!(bonds, i0)
        k += 16
    end
    @assert k==length(Js)
end

function push16!(bonds, i0)
    reds = i0 .+ (1:4)
    blacks = i0 .+ (5:8)
    for i=1:4
        for j=1:4
            push!(bonds, (reds[j], blacks[i]))
        end
    end
end

function chimera_bond!(bonds, Lx::Int, Ly::Int, J::AbstractVector)
    nj_red = (Ly-1)*4 + 16 * Ly
    println("Layer 1/$Lx")
    k = 0
    red_bond!(bonds, 1, Lx, Ly, J[k+1:k+nj_red])
    k += nj_red

    for j=2:Lx
        println("Layer $j/$Lx")
        # BLACK
        for iy=1:Ly
            i0 = index0(j-1, iy, Lx)
            j0 = index0(j, iy, Lx)
            for i = 1:4
                push!(bonds, (i0+4+i, j0+4+i))
                k += 1
            end
        end

        # Contract with RED
        rr = red_bond!(bonds, j, Lx, Ly, J[k+1:k+nj_red])
        k += nj_red
    end
    if length(J) !== k
        @warn "length of parameters is $(length(J)), used $k"
    end
    return bonds
end

if abspath(PROGRAM_FILE) == @__FILE__
    dump_JCs(Val(:randn))
    dump_Js(Val(:randn))
end

function i2cind(i::Int, Lx::Int, Ly::Int)
    cellid = mod1(i, 8)
    color = cellid <=4 ? :red : :black
    blockid = (i-1) ÷ 8
    iy = blockid ÷ Lx + 1
    ix = blockid - (iy-1)*Lx + 1
    ix, iy, cellid, color
end

using Test
@testset "assign labels" begin
    Lx, Ly = 4, 4
    Js = load_JC(Lx,Ly,Val(:randpm))
    labels = chimera_bond!([], Lx, Ly, Js)
    @test length(labels) == length(Js)
    indices = [Base.Iterators.Flatten(labels)...]
    for i=1:Lx*Ly*8
        ix, iy, cellid, color = i2cind(i, Lx, Ly)
        if (ix == 1 || ix == Lx) && color == :black
            @test count(==(i), indices) == 5
            continue
        end
        if (iy == 1 || iy == Ly) && color == :red
            @test count(==(i), indices) == 5
            continue
        end
        @test count(==(i), indices) == 6
    end
end
