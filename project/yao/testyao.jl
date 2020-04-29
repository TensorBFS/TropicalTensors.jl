include("spinglass.jl")

function first_row(L::Int, jtype::Val)
    I2 = spinglass_vertex_tensor(2)
    I3 = spinglass_vertex_tensor(3)
    ball = I2
    for j in 2:L-1
        ball = ein"(ab,bc),cde->ade"(ball, spinglass_bond_tensor(_get_J(jtype)), I3)
        ball = reshape(ball,:,size(ball,3))
    end
    ball = ein"(ab,bc),cd->ad"(ball, spinglass_bond_tensor(_get_J(jtype)), I2)
    vec(ball)
end

function one_more_stack(L::Int, ball, jtype::Val)
    get_TroB() = spinglass_bond_tensor(_get_J(jtype))
    I4 = spinglass_vertex_tensor(4)
    I3 = spinglass_vertex_tensor(3)
    ball = reshape(ball,2,:)
    ball = ein"(ab,ae),cde->cdb"(ball, get_TroB(), I3)
    ball = reshape(ball,2,2,2,:)
    for j in 2:L-1
        ball = ein"abcd,(be,(hc,efgh))->afgd"(ball, get_TroB(), get_TroB(),I4)
        ball = reshape(ball, size(ball,1)*size(ball,2),2,2,:)
    end
    ball = ein"abcd,(be,(hc,efh))->afd"(ball, get_TroB(), get_TroB(), I3)
    vec(ball)
end

function last_row(L, ball, jtype)
    I2 = spinglass_vertex_tensor(2)
    I3 = spinglass_vertex_tensor(3)
    get_TroB() = spinglass_bond_tensor(_get_J(jtype))
    ball = reshape(ball, 2, :)
    ball = ein"(ab,da),cd->cb"(ball, get_TroB(), I2)
    ball = reshape(ball,2,2,:)
    for j in 2:L-1
        ball = ein"abc,(ad,(fb,def))->ec"(ball, get_TroB(), get_TroB(), I3)
        ball = reshape(ball,2,2,:)
    end
    ball = ein"abc,(ad,(fb,df))->c"(ball, get_TroB(), get_TroB(), I2)
    return ball
end


@testset "yao" begin
    L = 10
    jtype = Val(:ferro)
    # Yao gates
    reg = ArrayReg(ones(Tropical{Float64}, 1<<L))
    G2 = matblock(spinglass_bond_tensor(1.0))
    G4 = matblock(spinglass_g4_tensor(1.0))

    # first row
    row1 = first_row(L, jtype)
    for i=1:L-1
        reg |> put(L, (i,i+1)=>G4)
    end
    @test row1 ≈ statevec(reg)

    # one more stack
    row2 = one_more_stack(L, row1, jtype)

    for i=1:L
        reg |> put(L, i=>G2)
    end
    for i=1:L-1
        reg |> put(L, (i,i+1)=>G4)
    end
    @test row2 ≈ statevec(reg)

    # last row
    # one more stack
    row3 = last_row(L, row2, jtype)
    for i=1:L
        reg |> put(L, i=>G2)
    end
    for i=1:L-1
        reg |> put(L, (i,i+1)=>G4)
    end
    @test row3[] ≈ sum(statevec(reg))
    @test spinglass_yao(10, Val(:ferro)) ≈ Tropical(180.0)
    @test spinglass_yao(10, Val(:ferro); usecuda=true) ≈ Tropical(180.0)
end
