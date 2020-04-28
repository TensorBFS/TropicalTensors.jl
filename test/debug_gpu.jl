using OMEinsum
using CuArrays
using TropicalTensors
CuArrays.allowscalar(false)

a = CuArray(Tropical.(ones(Float16,2,2,2^30,2)))
println("a generated")
b = CuArray(Tropical.(ones(Float16,2,2)))
#a = ein"abc,cd->abd"(a,b)
a[1,1,:,:] = ein"bc,cd->bd"(a[1,1,:,:],b)
println("stage 1")
a[1,2,:,:] = ein"bc,cd->bd"(a[1,2,:,:],b)
println("stage 2")
a[2,1,:,:] = ein"bc,cd->bd"(a[2,1,:,:],b)
println("stage 3")
a[2,2,:,:] = ein"bc,cd->bd"(a[2,2,:,:],b)
println("stage 4")




