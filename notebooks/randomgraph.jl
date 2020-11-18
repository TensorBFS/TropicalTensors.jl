### A Pluto.jl notebook ###
# v0.12.10

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
	Compose.set_default_graphic_size(14cm, 14cm)
	
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
		locs_x, locs_y = spring_layout(g)
		metas = [TensorMeta((0.5+0.5*locs_x[i], 0.5+0.5*locs_y[i]), string('A'+(i-1))) for i=1:n]
		TensorNetwork(tensors; metas=metas)
	end
	
	function spring_layout(g::AbstractGraph,
						   locs_x=2*rand(nv(g)).-1.0,
						   locs_y=2*rand(nv(g)).-1.0;
						   C=2.0,
						   MAXITER=100,
						   INITTEMP=2.0)

		nvg = nv(g)
		adj_matrix = adjacency_matrix(g)

		# The optimal distance bewteen vertices
		k = C * sqrt(4.0 / nvg)
		k² = k * k

		# Store forces and apply at end of iteration all at once
		force_x = zeros(nvg)
		force_y = zeros(nvg)

		# Iterate MAXITER times
		@inbounds for iter = 1:MAXITER
			# Calculate forces
			for i = 1:nvg
				force_vec_x = 0.0
				force_vec_y = 0.0
				for j = 1:nvg
					i == j && continue
					d_x = locs_x[j] - locs_x[i]
					d_y = locs_y[j] - locs_y[i]
					dist²  = (d_x * d_x) + (d_y * d_y)
					dist = sqrt(dist²)

					if !( iszero(adj_matrix[i,j]) && iszero(adj_matrix[j,i]) )
						# Attractive + repulsive force
						# F_d = dist² / k - k² / dist # original FR algorithm
						F_d = dist / k - k² / dist²
					else
						# Just repulsive
						# F_d = -k² / dist  # original FR algorithm
						F_d = -k² / dist²
					end
					force_vec_x += F_d*d_x
					force_vec_y += F_d*d_y
				end
				force_x[i] = force_vec_x
				force_y[i] = force_vec_y
			end
			# Cool down
			temp = INITTEMP / iter
			# Now apply them, but limit to temperature
			for i = 1:nvg
				fx = force_x[i]
				fy = force_y[i]
				force_mag  = sqrt((fx * fx) + (fy * fy))
				scale      = min(force_mag, temp) / force_mag
				locs_x[i] += force_x[i] * scale
				locs_y[i] += force_y[i] * scale
			end
		end

		# Scale to unit square
		min_x, max_x = minimum(locs_x), maximum(locs_x)
		min_y, max_y = minimum(locs_y), maximum(locs_y)
		function scaler(z, a, b)
			2.0*((z - a)/(b - a)) - 1.0
		end
		map!(z -> scaler(z, min_x, max_x), locs_x, locs_x)
		map!(z -> scaler(z, min_y, max_y), locs_y, locs_y)

		return locs_x, locs_y
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
	g = LightGraphs.SimpleGraph(length(neighbors) ÷ 2, [neighbors[:,i].+1 for i=1:n])
	locs_x, locs_y = spring_layout(g)
	metas = [TensorMeta((0.5+locs_x[i]/2, 0.5+locs_y[i]/2), string(i)) for i=1:n]
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
# ╠═3c7bb01c-2235-11eb-2c39-8bbcd2b9a08c
# ╠═c521ec84-2232-11eb-051e-bbbbfaf9c128
# ╠═84b2cb76-2240-11eb-0bf3-99e314c513be
# ╠═bcefc45c-2232-11eb-2615-21c5314bf0f3
# ╠═afec0844-2230-11eb-23ea-453a414ba792
# ╟─6b01ff60-217a-11eb-2a5c-13eb2716c727
# ╠═2d5dd844-2136-11eb-176f-e53e2ac97d82
# ╟─1aabc318-2137-11eb-33b7-cfcec36e631f
# ╟─72d7474c-2137-11eb-030b-eb894e4b22dc
# ╠═09ea2b64-2137-11eb-187a-fdf0f5b882ca
