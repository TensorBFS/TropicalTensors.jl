using CUDA, CuYao
using DelimitedFiles
device!(parse(Int, ARGS[1]))

include("dataload.jl")

function load_and_contract(::Type{T}, n::Int; seed::Int, usecuda=false, datafile="ising.hdf5") where T
	loadeddata = HDF5.h5open(TropicalTensors.project_relative_path("data", datafile), "r")
    instance = loadeddata["n$n"]["seed$seed"]
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
    tn = TensorNetwork(tensors)
    if usecuda
        tn = TropicalTensors.togpu(tn)
    end
    tree = build_tree(order)
    Array(TropicalTensors.contract(tn, tree).array)[]
end

function run(::Type{T}, n::Int; dataset) where T
    saveto = joinpath(@__DIR__, "$(dataset)_n$(n)_elsl.dat")
    elsl = zeros(T, 3, 100)
    datafile = dataset*".hdf5"
	loadeddata = HDF5.h5open(TropicalTensors.project_relative_path("data", datafile), "r")
    println("number of seeds = $(length(read(loadeddata["n$n"])))")
    t = @elapsed res = load_and_contract(T, n; seed=1, usecuda=true, datafile=dataset*".hdf5")
    for seed = 1:100
        try
            t = @elapsed res = load_and_contract(T, n; seed=seed, usecuda=true, datafile=datafile)
            println("seed = $seed, Δt = $t, maximum energy = $(res.n), degeneracy = $(res.c)")
            elsl[1,seed] = res.n/n
            elsl[2,seed] = log(res.c)/n
            elsl[3,seed] = t
        catch e
            println("Fail on seed $seed. $e")
        end
    end
    writedlm(saveto, elsl)
    return elsl
end

const n = parse(Int, ARGS[3])
@time run(Float32, n; dataset=ARGS[2])
