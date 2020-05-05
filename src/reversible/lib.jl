export igemm!, igemv!
export muleq, muleq_mul, muleq_add, tropical_muladd

export maxloc

# out initialized to one(Tropical)
@i function igemm!(out!::AbstractMatrix{T}, x::AbstractMatrix{T}, y::AbstractMatrix{T}) where T<:Tropical
	@safe size(x, 2) == size(y, 1) || throw(DimensionMismatch())
	@invcheckoff branch_keeper ← zeros(Bool, size(x,2))
	@invcheckoff for j=1:size(y,2)
		for i=1:size(x,1)
			el ← zero(T)
			@routine @inbounds for k=1:size(x,2)
				muleq(x[i,k], y[k,j])
				if (el.n < x[i,k].n, branch_keeper[k])
					FLIP(branch_keeper[k])
					NiLang.SWAP(el, x[i,k])
				end
			end
			@inbounds muleq(out![i,j], el)
			~@routine
		end
	end
	@invcheckoff branch_keeper → zeros(Bool, size(x,2))
end

@i function igemv!(out!::AbstractVector{T}, x::AbstractMatrix{T}, y::AbstractVector{T}, branch_keeper::AbstractVector{Bool}) where T<:Tropical
	@safe size(x, 2) == size(y, 1) || throw(DimensionMismatch())
	@invcheckoff for i=1:size(x,1)
		el ← zero(T)
		@routine @inbounds for k=1:size(x,2)
			muleq(x[i,k], y[k])
			if (el.n < x[i,k].n, branch_keeper[k])
				FLIP(branch_keeper[k])
				NiLang.SWAP(el, x[i,k])
			end
		end
		@inbounds muleq(out![i], el)
		~@routine
	end
end

NiLang.SWAP(a::Tropical{T}, b::Tropical{T}) where T = b, a
NiLang.SWAP(a::Array, b::Array) where T = b, a
# branch should be initialized to false.
@i @inline function tropical_muladd(out!::Tropical, x::Tropical, y::Tropical, branch)
	muleq(x, y)
	@invcheckoff if (out!.n < x.n, branch)
		FLIP(branch)
		NiLang.SWAP(out!, x)
	end
end

@i @inline function muleq_add(z::Tropical, x::Tropical, y::Tropical)
    @invcheckoff if (x.n > y.n, ~)
        z.n += identity(x.n)
    else
        z.n += identity(y.n)
    end
end

@i @inline function muleq(x::Tropical, y::Tropical)
    x.n += identity(y.n)
end

@i @inline function muleq_mul(out!::Tropical, x::Tropical, y::Tropical)
    x.n += x.n + y.n
end

maxloc(v::AbstractVector) = findmax(v)[2]

# AD wrappers
(_::Type{Inv{Tropical}})(x::Tropical) = x.n
NiLang.AD.GVar(x::Tropical) = Tropical(GVar(x.n))
NiLang.AD.grad(x::Tropical{<:GVar}) = Tropical(grad(x.n))
(_::Type{Inv{GVar}})(x::Tropical{<:GVar}) = Tropical((~GVar)(x.n))
