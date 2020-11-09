### A Pluto.jl notebook ###
# v0.12.7

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
	using PlutoMustache
	using TropicalTensors
	using Random
	using LightGraphs
	using Compose
	using HDF5
	Compose.set_default_graphic_size(10cm, 10cm)
	
	function rand_3regular_tn(::Type{T}, n; D=2) where T
		g = LightGraphs.random_regular_graph(n, 3)
		labels = 1:ne(g)
		arrays = [rand(T, fill(2, 3)...) for i=1:n]
		labels = [Int[] for i=1:n]
		for (k, e) in enumerate(edges(g))
			push!(labels[e.src], k)
			push!(labels[e.dst], k)
		end
		tensors = LabeledTensor.(arrays, labels)
		metas = [TensorMeta((0.5+0.5*cos(i/n*2π), 0.5+0.5*sin(i/n*2π)), string('A'+(i-1))) for i=1:n]
		TensorNetwork(tensors; metas=metas)
	end
end

# ╔═╡ cb989180-2216-11eb-2801-9fdfc4cec1d7
loadeddata = let
	f = HDF5.h5open(TropicalTensors.project_relative_path("data", "ising.hdf5"), "r")
	read(f)
end;

# ╔═╡ 3ad9d536-2217-11eb-0be0-4749c96f6ef5
md"""
instance size = $(@bind instance_size Select(collect(keys(loadeddata))))
"""

# ╔═╡ a17d9702-2217-11eb-2213-a15da2f07bed
loadeddata_sized = loadeddata[instance_size];

# ╔═╡ 8f652f7e-2217-11eb-38aa-73b178bef9cd
md"""
seed = $(@bind instance_seed Select(collect(keys(loadeddata_sized))))
"""

# ╔═╡ 2316240e-2230-11eb-0b78-4182b35b1b77
loadeddata_sized[instance_seed] |> keys

# ╔═╡ 2717c314-2217-11eb-03b4-b3ee8f221a02
begin
	instance = loadeddata_sized[instance_seed]
	arrays = instance["tensors"]
	entros = instance["entros"]
	order = instance["order"]
	neighbors = instance["labels"]
end;

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

# ╔═╡ 3c7bb01c-2235-11eb-2c39-8bbcd2b9a08c
tn = let
	n = size(arrays, 4)
	labels = tensor_labels(neighbors)
	tensors = map(1:n) do i
		arr = Float64.(permutedims(arrays[:,:,:,i], (3,2,1)))
		degen = round.(exp.(permutedims(entros[:,:,:,i], (3,2,1))))
		LabeledTensor(CountingTropical{Float64}.(arr, degen), labels[:,i]) 
	end
	# circle layout
	metas = [TensorMeta((0.5+0.5*cos(i/n*2π), 0.5+0.5*sin(i/n*2π)), string(i)) for i=1:n]
	TensorNetwork(tensors, metas=metas)
end

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

# ╔═╡ 84b2cb76-2240-11eb-0bf3-99e314c513be
neighbors

# ╔═╡ bcefc45c-2232-11eb-2615-21c5314bf0f3
tree = build_tree(order)

# ╔═╡ afec0844-2230-11eb-23ea-453a414ba792
TropicalTensors.contract(tn, tree).array[]

# ╔═╡ 6b01ff60-217a-11eb-2a5c-13eb2716c727
md"number of nodes = $(@bind nnodes Slider(2:2:40, default=16))"

# ╔═╡ 2d5dd844-2136-11eb-176f-e53e2ac97d82
tnr = rand_3regular_tn(Float64, nnodes; D=2)

# ╔═╡ 1aabc318-2137-11eb-33b7-cfcec36e631f
@bind step Button(label="contract")

# ╔═╡ 72d7474c-2137-11eb-030b-eb894e4b22dc
@regstate step

# ╔═╡ 09ea2b64-2137-11eb-187a-fdf0f5b882ca
@when step begin
	length(tn.tensors) > 1 && contract!(tn, rand(vcat([t.labels for t in tn.tensors]...)))
	tn
end

# ╔═╡ Cell order:
# ╠═c7118f3c-209e-11eb-1f79-bdc2dd17aea0
# ╟─cb989180-2216-11eb-2801-9fdfc4cec1d7
# ╟─3ad9d536-2217-11eb-0be0-4749c96f6ef5
# ╠═a17d9702-2217-11eb-2213-a15da2f07bed
# ╟─8f652f7e-2217-11eb-38aa-73b178bef9cd
# ╠═2316240e-2230-11eb-0b78-4182b35b1b77
# ╠═2717c314-2217-11eb-03b4-b3ee8f221a02
# ╠═3d04dd42-223a-11eb-14f7-9d3d29d35349
# ╟─3c7bb01c-2235-11eb-2c39-8bbcd2b9a08c
# ╠═c521ec84-2232-11eb-051e-bbbbfaf9c128
# ╠═84b2cb76-2240-11eb-0bf3-99e314c513be
# ╠═bcefc45c-2232-11eb-2615-21c5314bf0f3
# ╠═afec0844-2230-11eb-23ea-453a414ba792
# ╟─6b01ff60-217a-11eb-2a5c-13eb2716c727
# ╠═2d5dd844-2136-11eb-176f-e53e2ac97d82
# ╟─1aabc318-2137-11eb-33b7-cfcec36e631f
# ╟─72d7474c-2137-11eb-030b-eb894e4b22dc
# ╠═09ea2b64-2137-11eb-187a-fdf0f5b882ca
