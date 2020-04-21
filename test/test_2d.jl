using EliminateGraphs
using TropicalTensors
using TropicalTensors: inferier_table
using OMEinsum

function get_TroI()
    I2 = zeros(2,2)
    I2 .= -Inf
    I2[1,1] = 0.0
    I2[2,2] = 0.0

    I3 = zeros(2,2,2)
    I3 .= -Inf
    I3[1,1,1] = 0.0
    I3[2,2,2] = 0.0

    I4 = zeros(2,2,2,2)
    I4 .= -Inf
    I4[1,1,1,1] = 0.0
    I4[2,2,2,2] = 0.0

    return Tropical.(I2),Tropical.(I3),Tropical.(I4)
end

function get_TroB(Jtype::String)
    if (Jtype == "ferro")
        Jij = 1.0
    elseif (Jtype == "randn")
        Jij = randn()
    elseif (Jtype == "rand")
        Jij = rand()
    else
        println("wrong Jij type")
        exit()
    end
    return Tropical.([Jij -Jij;-Jij Jij])
end

function gen_mine_2d(L::Int, Jtype::String)
    I2,I3,I4 = get_TroI()
    """The first row"""
    ball = I2
    for j in 2:L-1
        ball = ein"ab,bc,cde->ade"(ball,get_TroB(Jtype),I3)
        ball = reshape( ball,:,size(ball,3) )
    end
    ball = ein"ab,bc,cd->ad"(ball,get_TroB(Jtype),I2)

    """ From the second row to the second row from bottom """
    for i in 2:L-1
        ball = reshape(ball,2,:)
        ball = ein"ab,ae,cde->cdb"(ball,get_TroB(Jtype),I3)
        ball = reshape(ball,2,2,2,:)
        for j in 2:L-1
            ball = ein"abcd,be,hc,efgh->afgd"(ball,get_TroB(Jtype),get_TroB(Jtype),I4)
            ball = reshape(ball, size(ball,1)*size(ball,2),2,2,:)
        end
        ball = ein"abcd,be,hc,efh->afd"(ball,get_TroB(Jtype),get_TroB(Jtype),I3)
        ball = reshape(ball,2,:)
    end

    """ the last row """
    ball = ein"ab,da,cd->cb"(ball,get_TroB(Jtype),I2)
    ball = reshape(ball,2,2,:)
    for j in 2:L-1
        ball = ein"abc,ad,fb,def->ec"(ball,get_TroB(Jtype),get_TroB(Jtype),I3)
        ball = reshape(ball,2,2,:)
    end
    ball = ein"abc,ad,fb,df->c"(ball,get_TroB(Jtype),get_TroB(Jtype),I2)
    return ball
end

Jtype = "ferro"
Jtype = "randn"
L=6
e = gen_mine_2d(L,Jtype)
println("mine=",e)
