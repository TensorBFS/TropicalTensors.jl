using CUDA
device!(parse(Int, get(ARGS, 1, "0")))

using TropicalTensors
using HDF5

function tensor_labels(neighbors)
	k=0
	nbs = copy(neighbors)
	labels = zeros(Int, size(nbs)...)
	for j=1:size(nbs, 2)
		for i=1:size(nbs, 1)
			j2 = nbs[i, j] + 1
			if j2 != 0
				k += 1
				labels[i, j] = k
				nbs[i, j] = -1
				i2 = indexin(j-1, nbs[:,j2])[]
				labels[i2, j2] = k
				nbs[i2, j2] = -1
			end
		end
	end
	@assert count(labels .== 0) == 0 && length(unique(labels)) == length(labels)÷2
	labels
end

function build_tree(order)
	N = size(order, 2)+1
	ids = collect(Any, 1:N)
	for i=1:N-1
		ta, tb = order[:,i]
		ids[ta+1] = ContractionTree(ids[ta+1], ids[tb+1])
		ids[tb+1] = nothing
	end
	ids[1]
end

function panzhang(::Type{T}, n::Int; seed::Int, usecuda=false) where T
	loadeddata = HDF5.h5open(TropicalTensors.project_relative_path("data", "ising.hdf5"), "r")
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
		degen = round.(T, exp.(permutedims(entros[:,:,:,i], (3,2,1))))
		LabeledTensor(CountingTropical{T}.(arr, degen), labels[:,i]) 
	end
	# circle layout
	metas = [TensorMeta((0.5+0.5*cos(i/n*2π), 0.5+0.5*sin(i/n*2π)), string(i)) for i=1:n]
    tn = TensorNetwork(tensors, metas=metas)
    if usecuda
        tn = togpu(tn)
    end
    tree = build_tree(order)
    TropicalTensors.contract(tn, tree).array[]
end

@time panzhang(160; seed=97, usecuda=true)