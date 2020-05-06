export igemm!, igemv!, isum
export muleq, muleq_mul, muleq_add, tropical_muladd

export maxloc

const TropicalG{T,TG} = Tropical{GVar{T,TG}}

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

NiLang.SWAP(a, b) where T = b, a
# branch should be initialized to false.
@i @inline function tropical_muladd(out!::Tropical, x::Tropical, y::Tropical, branch)
	muleq(x, y)
	@invcheckoff if (out!.n < x.n, branch)
		FLIP(branch)
		NiLang.SWAP(out!, x)
	end
end

@i function isum(out!::T, v::AbstractArray{T}) where T<:Tropical
	@routine @invcheckoff begin
		branch_keeper ← zeros(Bool, length(v))
		anc ← zero(T)
		for i = 1:length(v)
			if (anc.n < v[i].n, branch_keeper[i])
				FLIP(branch_keeper[i])
				NiLang.SWAP(anc, v[i])
			end
		end
	end
	out!.n += identity(anc.n)
	~@routine
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
NiLang.AD.GVar(x::Tropical) = Tropical(GVar(x.n, zero(x.n)))
Base.convert(::Type{GVar{T,T}}, x::TX) where {T<:Number, TX<:Number} = GVar(T(x))
Base.convert(::Type{GVar{T,T}}, x::GVar{T,T}) where {T<:Number} = x
#(::Type{GVar{T,T}})(x::TX) where {T<:Number, TX<:Number} = GVar(T(x))
NiLang.AD.grad(x::TropicalG) = Tropical(grad(x.n))
(_::Type{Inv{GVar}})(x::TropicalG) = Tropical((~GVar)(x.n))
Base.isfinite(x::GVar) = isfinite(x.x) && isfinite(x.g)
function NiLang.loaddata(::Type{Array{TropicalG{T,GT},N}}, data::Array{Tropical{T},N}) where {T,GT,N}
    GVar.(data)
end
import NiLang.NiLangCore: deanc
#function deanc(x::Array{Tropical{GVar{T,T}},N}, val::Array{Tropical{T},N}) where {T,N}
#    deanc.(x, val)
#end
#function deanc(x::Tropical{GVar{T,T}}, val::Tropical{T}) where {T,N}
#    deanc(x.n, val.n)
#end

function deanc(x::T, val::T) where {T<:Array}
    deanc.(x, val)
end
function deanc(x::T, val::T) where T<:Tropical
    deanc(x.n, val.n)
end

Base.isapprox(x::GVar, y::GVar; kwargs...) = isapprox(x.x, y.x; kwargs...) && isapprox(x.g, y.g; kwargs...)

Base.one(x::XT) where XT<:TropicalG{T,GT} where {T,GT} = one(TropicalG{T,GT})
Base.one(::Type{TropicalG{T,GT}}) where {T,GT} = Tropical(GVar(zero(T), zero(GT)))
Base.zero(x::XT) where XT<:TropicalG{T,GT} where {T,GT} =zero(TropicalG{T,GT})
Base.zero(::Type{TropicalG{T,GT}}) where {T<:AbstractFloat,GT} = Tropical(GVar(typemin(T)/2, zero(GT)))
Base.zero(::Type{TropicalG{T,GT}}) where {T<:Integer,GT} = Tropical(GVar(typemin(T)÷2, zero(GT)))
