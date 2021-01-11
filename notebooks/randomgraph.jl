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

# ╔═╡ 6b01ff60-217a-11eb-2a5c-13eb2716c727
md"number of nodes = $(@bind nnodes Slider(2:2:40, default=16))"

# ╔═╡ 2d5dd844-2136-11eb-176f-e53e2ac97d82
tn = rand_3regular_tn(Float64, nnodes; D=2)

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
# ╟─6b01ff60-217a-11eb-2a5c-13eb2716c727
# ╠═2d5dd844-2136-11eb-176f-e53e2ac97d82
# ╟─1aabc318-2137-11eb-33b7-cfcec36e631f
# ╟─72d7474c-2137-11eb-030b-eb894e4b22dc
# ╠═09ea2b64-2137-11eb-187a-fdf0f5b882ca
