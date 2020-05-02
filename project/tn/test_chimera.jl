using TropicalTensors
using OMEinsum
using CuArrays
CuArrays.allowscalar(false)


function get_C(dtype)
    L=16
    return CuArray(Tropical.(ones(dtype,2^L,2^L)))
end

function test_cont()
    t0 = time()
    #dtype = Float32
    dtype = Float16

    print("generatig C1 and C2... ")
    C1 = get_C(dtype)
    C2 = get_C(dtype)
    println("done, ",time()-t0)
#    println("try vector")
#    C1 = C1[1,:]
#    println("shape ",size(C1)," ",size(C2))
#    C1 = ein"j,jk->k"(C1,C2)
#    print("done",time()-t0)

    chunks = 2^9
    C1  = reshape(C1,chunks,Int(size(C1,1)/chunks),:)
    for i in 1:chunks
        t1 = time()
        print("C1*C2 ", i,"/",chunks," ... ")
        C1[i,:,:] = ein"ij,jk->ik"(C1[i,:,:],C2)
        println("done",time()-t1, " Sec. total time",time()-t0," Sec.")
    end

#    print("C1*C2... ")
#    C1 = ein"ij,jk->ik"(C1,C2)
#    println("done, ",time()-t0)
#    

#    print("generatig C3... ")
#    C2 = get_C(dtype)
#    println("done, ",time()-t0)
#    print("C1*C3... ")
#    C1 = ein"ij,jk->ik"(C1,C2)
#    println("done, ",time()-t0)
#    print("generatig C4... ")
#    C2 = get_C(dtype)
#    print("overlap... ")
#    C1 = ein"ij,ji->"(C1,C2)
#    println("done, ",time()-t0)
end

test_cont()
