gens(L::SimpleNumField{T}) where {T} = [gen(L)]

function gen(L::SimpleNumField{T}, i::Int) where {T}
   i == 1 || error("index must be 1")
   return gen(L)
end

function Base.getindex(L::SimpleNumField{T}, i::Int) where {T}
   if i == 0
      return one(L)
   elseif i == 1
      return gen(L)
   else
      error("index has to be 0 or 1")
   end
end

number_of_generators(L::SimpleNumField{T}) where {T} = 1

is_unit(a::NumFieldElem) = !iszero(a)

canonical_unit(a::NumFieldElem) = a

characteristic(F::NumField) = 0

promote_rule(::Type{T}, ::Type{S}) where {S<:NumFieldElem,T<:Integer} = S

promote_rule(::Type{S}, ::Type{T}) where {S<:NumFieldElem,T<:Integer} = S
