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

# ╔═╡ c5eda0dc-215e-11eb-22d1-d76f1b9e047c
using Compose, Viznet, PlutoUI, TropicalTensors

# ╔═╡ 2824ef46-2165-11eb-23c9-a3335e6a6dba
md"## Define a lattice"

# ╔═╡ 47cbb81c-2169-11eb-0b07-9ff0d4a4c6b9
@bind lattice_type Select(["square", "chimera"]; default="square")

# ╔═╡ 9725485e-2161-11eb-2cf2-99d30b6a6dd6
md"""
``L_x`` = $(@bind Lx Slider(1:(lattice_type == "square" ? 22 : 6); default=(lattice_type == "square" ? 10 : 4), show_value=true))
"""

# ╔═╡ b320dc3a-2161-11eb-3edf-a1521d67bae0
md"""
``L_y`` = $(@bind Ly Slider(1:(lattice_type == "square" ? 22 : 6); default=(lattice_type == "square" ? 10 : 4), show_value=true))
"""

# ╔═╡ c9a88462-2161-11eb-186d-c9c32bbc4af5
if Lx > Ly
	Compose.set_default_graphic_size(14cm, 14cm/Lx*Ly)
else
	Compose.set_default_graphic_size(14cm/Ly*Lx, 14cm)
end

# ╔═╡ 126df144-215f-11eb-116e-a179ec34b675
lt = lattice_type == "square" ? SquareLattice(Lx, Ly) : ChimeraLattice(Lx, Ly)

# ╔═╡ 370cda48-2165-11eb-1514-d31e6fb4eb55
md"## Define the coupling"

# ╔═╡ b7cdbb3a-2162-11eb-0f41-6d89bcb699a6
md"""
coupling type = $(@bind coupling_type Select(["±1", "normal"]; default="±1"))
"""

# ╔═╡ 2cd5c642-2164-11eb-3e69-fd9f6da760fb
Js = if coupling_type == "±1"
	rand([-1.0,1], length(sgbonds(lt)))
elseif coupling_type == "normal"
	randn(length(sgbonds(lt)))
end;

# ╔═╡ 64f4a1ca-2163-11eb-20e8-cb6cfef00ec5
md"""
magnetic field = $(@bind magnetic_field Select(["abscent", "±1", "normal"]; default="abscent"))
"""

# ╔═╡ 1a3456ac-215f-11eb-3c67-39ef0ab19296
hs = if magnetic_field == "±1"
	rand([-1.0,1], length(sgvertices(lt)))
elseif magnetic_field == "normal"
	randn(length(sgvertices(lt)))
elseif magnetic_field == "abscent"
	zeros(Float64, length(sgvertices(lt)));
end;

# ╔═╡ 46644ab6-215f-11eb-08f5-378294e5fe24
# yellow -> negative coupling
# blue -> positive coupling
sg = Spinglass(lt, Js, hs)

# ╔═╡ 4c38a7f6-2165-11eb-15fd-7559fe9bbe86
md"## Find the optimal configuration"

# ╔═╡ 45112fb0-2166-11eb-3617-450f9d427a5e
res = solve(CountingTropical{Float64}, sg; usecuda=false)

# ╔═╡ 6ff289d6-2166-11eb-3a5c-4f1e686b2810
md"The minimum energy is $(res.n), degenercy is $(res.c)"

# ╔═╡ 3f5c4dba-2162-11eb-14b9-632cdafa348d
optimal_config = Reversible.opt_config(sg)

# ╔═╡ Cell order:
# ╠═c5eda0dc-215e-11eb-22d1-d76f1b9e047c
# ╟─2824ef46-2165-11eb-23c9-a3335e6a6dba
# ╟─47cbb81c-2169-11eb-0b07-9ff0d4a4c6b9
# ╟─9725485e-2161-11eb-2cf2-99d30b6a6dd6
# ╟─b320dc3a-2161-11eb-3edf-a1521d67bae0
# ╟─c9a88462-2161-11eb-186d-c9c32bbc4af5
# ╠═126df144-215f-11eb-116e-a179ec34b675
# ╟─370cda48-2165-11eb-1514-d31e6fb4eb55
# ╟─b7cdbb3a-2162-11eb-0f41-6d89bcb699a6
# ╠═2cd5c642-2164-11eb-3e69-fd9f6da760fb
# ╟─64f4a1ca-2163-11eb-20e8-cb6cfef00ec5
# ╠═1a3456ac-215f-11eb-3c67-39ef0ab19296
# ╠═46644ab6-215f-11eb-08f5-378294e5fe24
# ╟─4c38a7f6-2165-11eb-15fd-7559fe9bbe86
# ╠═45112fb0-2166-11eb-3617-450f9d427a5e
# ╟─6ff289d6-2166-11eb-3a5c-4f1e686b2810
# ╠═3f5c4dba-2162-11eb-14b9-632cdafa348d