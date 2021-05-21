### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ c7118f3c-209e-11eb-1f79-bdc2dd17aea0
begin
	using Revise
	using PlutoUI
	using TropicalTensors
	using Random
	using LightGraphs, Compose
	using SimpleTensorNetworks
	using HDF5
	Compose.set_default_graphic_size(14cm, 14cm)
end

# ╔═╡ eccb04de-5441-11eb-3825-fde87c880d75
md"# Ising sipnglass and 2-SAT counting on a 3 regular graph"

# ╔═╡ d44aedc0-5441-11eb-2aad-617e9def7361
md"## Load the graph instances from the disk"

# ╔═╡ d56772da-5448-11eb-24f9-2b00feecb187
md"First, we load both the graph instance and the contraction order the from the disk. The contraction order project in Julia is still work in progress, so we are generating the data with python. For someone who are interested in contraction orders, please file an issue in our [github repo](https://github.com/TensorBFS/TropicalTensors.jl)."

# ╔═╡ cb989180-2216-11eb-2801-9fdfc4cec1d7
loadeddata = let
	f = HDF5.h5open(TropicalTensors.project_relative_path("data", "ising.hdf5"), "r")
	read(f)
end;

# ╔═╡ 3ad9d536-2217-11eb-0be0-4749c96f6ef5
md"""
instance size = $(@bind instance_size Select(sort!(collect(keys(loadeddata)), by=x->parse(Int, x[2:end]))))
"""

# ╔═╡ a17d9702-2217-11eb-2213-a15da2f07bed
loadeddata_sized = loadeddata[instance_size];

# ╔═╡ 8f652f7e-2217-11eb-38aa-73b178bef9cd
md"""
seed = $(@bind instance_seed Select(sort!(collect(keys(loadeddata_sized)), by=x->parse(Int,x[5:end]))))
"""

# ╔═╡ 5a2ac51a-5446-11eb-39ba-f58535fe234a
md"Each loaded dict records information about a tensor network"

# ╔═╡ 2316240e-2230-11eb-0b78-4182b35b1b77
loadeddata_sized[instance_seed] |> keys

# ╔═╡ d380b3b6-5446-11eb-2090-1fe2eaf046c8
md"In the following, we want to construct a tensor network form these loaded data"

# ╔═╡ 2717c314-2217-11eb-03b4-b3ee8f221a02
begin
	instance = loadeddata_sized[instance_seed]
	arrays = instance["tensors"]
	entros = instance["entros"]
	order = instance["order"]
	neighbors = instance["labels"]   # these are neighbors, not labels
end;

# ╔═╡ 9391674c-5445-11eb-2ae9-355ce96941b1
md"Define the function to assign labels to each tensor, the original data records the neighbors rather than tensor lebels."

# ╔═╡ 3d04dd42-223a-11eb-14f7-9d3d29d35349
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

# ╔═╡ a60db2c2-5445-11eb-0f22-b5a145883689
md"Build the tensor network, and permute the tensor dimensions (the original array data is generated from python in C order)."

# ╔═╡ 3c7bb01c-2235-11eb-2c39-8bbcd2b9a08c
tn = let
	n = size(arrays, 4)
	labels = tensor_labels(neighbors)
	tensors = map(1:n) do i
		# from C order to Julia's column major order
		arr = Float64.(permutedims(arrays[:,:,:,i], (3,2,1)))
		# from entropy to 
		degen = round.(exp.(permutedims(entros[:,:,:,i], (3,2,1))))
		# the first parameter is the tensor element type, where `CountingTropical{Float64}` is a tropical number type that do the counting while computing. The storage is 64 bit floating point numbers. One can also use `Tropical{Float64}` for computing ground state. The reason for not using `Int` is because the degeneracy can be very huge.
		LabeledTensor(CountingTropicalF64.(arr, degen), labels[:,i]) 
	end
	TensorNetwork(tensors)
end

# ╔═╡ 17092d36-5445-11eb-3952-9973673116b4
md"## Build the tree for the contraction order and do the contraction"

# ╔═╡ 036130fa-5448-11eb-32cb-41eb39edc4ad
md"The contraction order is recorded in a different way as the tree representation as we need for contraction, so we need to do a transformation."

# ╔═╡ c521ec84-2232-11eb-051e-bbbbfaf9c128
function build_tree(order)
	N = length(tn.tensors)
	ids = collect(Any, 1:N)
	for i=1:N-1
		ta, tb = order[:,i]
		ids[ta+1] = ContractionTree(ids[ta+1], ids[tb+1])
		ids[tb+1] = nothing
	end
	ids[1]
end

# ╔═╡ bcefc45c-2232-11eb-2615-21c5314bf0f3
tree = build_tree(order)

# ╔═╡ afec0844-2230-11eb-23ea-453a414ba792
res = TropicalTensors.contract(tn, tree).array[]

# ╔═╡ ee2a0b00-5445-11eb-0146-75121bef8e47
md"The optimal is $(res.n), with degeneracy is $(res.c)"

# ╔═╡ Cell order:
# ╟─eccb04de-5441-11eb-3825-fde87c880d75
# ╠═c7118f3c-209e-11eb-1f79-bdc2dd17aea0
# ╟─d44aedc0-5441-11eb-2aad-617e9def7361
# ╟─d56772da-5448-11eb-24f9-2b00feecb187
# ╠═cb989180-2216-11eb-2801-9fdfc4cec1d7
# ╟─3ad9d536-2217-11eb-0be0-4749c96f6ef5
# ╠═a17d9702-2217-11eb-2213-a15da2f07bed
# ╟─8f652f7e-2217-11eb-38aa-73b178bef9cd
# ╟─5a2ac51a-5446-11eb-39ba-f58535fe234a
# ╠═2316240e-2230-11eb-0b78-4182b35b1b77
# ╟─d380b3b6-5446-11eb-2090-1fe2eaf046c8
# ╠═2717c314-2217-11eb-03b4-b3ee8f221a02
# ╟─9391674c-5445-11eb-2ae9-355ce96941b1
# ╠═3d04dd42-223a-11eb-14f7-9d3d29d35349
# ╟─a60db2c2-5445-11eb-0f22-b5a145883689
# ╠═3c7bb01c-2235-11eb-2c39-8bbcd2b9a08c
# ╟─17092d36-5445-11eb-3952-9973673116b4
# ╟─036130fa-5448-11eb-32cb-41eb39edc4ad
# ╠═c521ec84-2232-11eb-051e-bbbbfaf9c128
# ╠═bcefc45c-2232-11eb-2615-21c5314bf0f3
# ╠═afec0844-2230-11eb-23ea-453a414ba792
# ╟─ee2a0b00-5445-11eb-0146-75121bef8e47
