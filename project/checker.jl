using Test
include("dataload.jl")

function panzhang_check(::Type{T}, n::Int; seed::Int, datafile="ising_test.hdf5") where T
	loadeddata = HDF5.h5open(TropicalTensors.project_relative_path("data", datafile), "r")
    instance = loadeddata["n$n"]["seed$seed"]
    #return instance
    arrays = read(instance, "tensors")
    entros = read(instance, "entros")
    order = read(instance, "order")
    neighbors = read(instance, "labels")
    E = read(instance, "E")
    S = read(instance, "S")
    close(loadeddata)

    # build tensor network
	labels = tensor_labels(neighbors)
	tensors = map(1:n) do i
		arr = T.(permutedims(arrays[:,:,:,i], (3,2,1)))
		degen = T.(round.(exp.(permutedims(entros[:,:,:,i], (3,2,1)))))
		LabeledTensor(CountingTropical{T}.(arr, degen), labels[:,i]) 
	end
	# circle layout
	metas = [TensorMeta((0.5+0.5*cos(i/n*2π), 0.5+0.5*sin(i/n*2π)), string(i)) for i=1:n]
    tn = TensorNetwork(tensors, metas=metas)
    tree = build_tree(order)
    res = Array(TropicalTensors.contract(tn, tree).array)[]
    println("Mine: E = $(res.n/n), S = $(log(res.c)/n)")
    println("PanZhang: E = $E, S = $S")
    @test isapprox(res.n, E * n; atol=1e-2)
    @test isapprox(res.c, exp(S*n), atol=1e-2)
end


@testset "checking" begin
    for n=[10,16]
        for seed = 1:3
            panzhang_check(Float64, n; seed=seed, datafile="ising_test_n1016_seed123.hdf5")
            panzhang_check(Float64, n; seed=seed, datafile="2sat_test_n1016_seed123.hdf5")
        end
    end
end
