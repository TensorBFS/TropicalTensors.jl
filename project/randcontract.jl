using CUDA, CuYao
device!(parse(Int, ARGS[1]))

include("panzhangreader.jl")

function panzhang(::Type{T}, n::Int; seed::Int, usecuda=false, datafile="ising.hdf5") where T
	loadeddata = HDF5.h5open(TropicalTensors.project_relative_path("data", datafile), "r")
    instance = loadeddata["n$n"]["seed$seed"]
    #return instance
    arrays = read(instance, "tensors")
    entros = read(instance, "entros")
    order = read(instance, "order")
    neighbors = read(instance, "labels")
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
    if usecuda
        tn = TropicalTensors.togpu(tn)
    end
    tree = build_tree(order)
    Array(TropicalTensors.contract(tn, tree).array)[]
end

function run(n::Int; saveto)
    elsl = zeros(2, 100)
    for seed = 1:100
        res = @time panzhang(Float64, n; seed=97, usecuda=true)
        @show res
        elsl[1,seed] = res.n/n
        elsl[2,seed] = log(res.c)/n
    end
    return elsl
end

const n = 100
run(n; saveto=joinpath(@__DIR__, "n$(n)_elsl.dat"))
