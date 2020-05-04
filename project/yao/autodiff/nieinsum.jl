# NOTE: this is slow!

"""
A naive implementation of `einsum!`
    * `ixs`: input tensor indices,
    * `xs`: input tensors,
    * `iy`: output tensor indices,
    * `y`: accumulated tensor, notice it is initialized to 0 as output!
"""
@i function naive_einsum!(code::EinCode{ixs, iy}, xs, y::AbstractArray{T}) where {ixs, iy, NO,T<:Tropical}
	@routine @invcheckoff begin
	    # outer legs and inner legs
	    outer_indices ← unique(iy)
	    inner_indices ← setdiff(TupleTools.vcat(ixs...), outer_indices)

	    # find size for each leg
	    all_indices ← TupleTools.vcat(ixs..., iy)
	    all_sizes ← TupleTools.vcat(size.(xs)..., size(y))
	    outer_sizes ← [map(i->all_sizes[i], indexin(outer_indices, [all_indices...]))...]
	    inner_sizes ← [map(i->all_sizes[i], indexin(inner_indices, [all_indices...]))...]

	    # cartesian indices for outer and inner legs
	    outer_ci ← CartesianIndices((outer_sizes...,))
	    inner_ci ← CartesianIndices((inner_sizes...,))

	    # for indexing tensors (leg binding)
	    indices ← (outer_indices..., inner_indices...)
	    locs_xs ← map(ix->map(i->findfirst(isequal(i), indices), ix), ixs)
    	locs_y ← map(i->findfirst(isequal(i), outer_indices), iy)
	end
    loop!(locs_xs, xs, locs_y, y, outer_ci, inner_ci)
	~@routine
end

"""take an index subset from `ind`"""
index_map(ind::CartesianIndex, locs::Tuple) = CartesianIndex(TupleTools.getindices(Tuple(ind), locs))

"""
loop and accumulate products to y, the GPU version, the CPU version.
"""
@i function loop!(locs_xs::NTuple{N}, xs::NTuple{N, AbstractArray}, locs_y, y::AbstractArray{T}, outer_ci::CartesianIndices, inner_ci::CartesianIndices) where {N, T<:Tropical}
    @invcheckoff @inbounds for i in outer_ci
        ind_y ← outer_ci[i]
        iy ← index_map(ind_y, locs_y)
		branch_keeper ← zeros(Bool, length(inner_ci))
		pl ← ones(T, length(inner_ci))
		el ← zero(T)
		@routine for ind_x in inner_ci
            ind_xy ← CartesianIndex(TupleTools.vcat(ind_y.I, ind_x.I))
            #y[iy] += map_prod(T, xs, ind_xy, locs_xs)
            for I=1:N
                muleq(pl[ind_x], xs[I][index_map(ind_xy, locs_xs[I])])
            end
			if (el.n < pl[ind_x].n, branch_keeper[ind_x])
				FLIP(branch_keeper[ind_x])
				NiLang.SWAP(el, pl[ind_x])
			end
        end
		muleq(y[iy], el)
		~@routine
    end
end

using Test
N = 3
a = Tropical.(randn(N, N))
b = Tropical.(randn(N, N))
c = Tropical.(zeros(N, N))
naive_einsum!(ein"ij,jk->ik", (a, b), c)
@test c ≈ a*b
