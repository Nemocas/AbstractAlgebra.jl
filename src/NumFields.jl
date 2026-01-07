gens(L::SimpleNumField{T}) where {T} = [gen(L)]

function gen(L::SimpleNumField{T}, i::Int) where {T}
   i == 1 || error("index must be 1")
   return gen(L)
end

number_of_generators(L::SimpleNumField{T}) where {T} = 1

characteristic(F::NumField) = 0

promote_rule(::Type{T}, ::Type{S}) where {S<:NumFieldElem,T<:Integer} = S

promote_rule(::Type{S}, ::Type{T}) where {S<:NumFieldElem,T<:Integer} = S
