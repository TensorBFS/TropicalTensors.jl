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
	using PlutoMustache
	using TropicalTensors
	using Random
	using LightGraphs
	using Compose, Viznet
	using SimpleTensorNetworks
	using HDF5
	Compose.set_default_graphic_size(14cm, 14cm)
end

# ╔═╡ 73cff25e-5434-11eb-0147-21ac9dbcb8b5
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
	TensorNetwork(tensors)
end

# ╔═╡ 6b01ff60-217a-11eb-2a5c-13eb2716c727
md"number of nodes = $(@bind nnodes Slider(2:2:40, default=16))"

# ╔═╡ 2d5dd844-2136-11eb-176f-e53e2ac97d82
tnr = rand_3regular_tn(Float64, nnodes; D=2)

# ╔═╡ 8d7cb70c-5441-11eb-35cb-e551b9166fd6
md"contract the graph by a random order"

# ╔═╡ 1aabc318-2137-11eb-33b7-cfcec36e631f
@bind step Button(label="step")

# ╔═╡ 72d7474c-2137-11eb-030b-eb894e4b22dc
@regstate step

# ╔═╡ 09ea2b64-2137-11eb-187a-fdf0f5b882ca
@when step begin
	length(tnr.tensors) > 1 && SimpleTensorNetworks.contract_label!(tnr, rand(vcat([t.labels for t in tnr.tensors]...)))
	tnr
end

# ╔═╡ Cell order:
# ╠═c7118f3c-209e-11eb-1f79-bdc2dd17aea0
# ╠═73cff25e-5434-11eb-0147-21ac9dbcb8b5
# ╟─6b01ff60-217a-11eb-2a5c-13eb2716c727
# ╠═2d5dd844-2136-11eb-176f-e53e2ac97d82
# ╟─8d7cb70c-5441-11eb-35cb-e551b9166fd6
# ╟─1aabc318-2137-11eb-33b7-cfcec36e631f
# ╟─72d7474c-2137-11eb-030b-eb894e4b22dc
# ╠═09ea2b64-2137-11eb-187a-fdf0f5b882ca
