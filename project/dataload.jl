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