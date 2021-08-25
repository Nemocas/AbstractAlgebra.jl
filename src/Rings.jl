promote_rule(::Type{T}, ::Type{T}) where T <: RingElement = T

function promote_rule_sym(::Type{T}, ::Type{S}) where {T, S}
   U = promote_rule(T, S)
   if U !== Union{}
      return U
   else
      UU = promote_rule(S, T)
      return UU
   end
end

@inline function try_promote(x::S, y::T) where {S <: RingElem, T <: RingElem}
   U = promote_rule_sym(S, T)
   if S === U
      return true, x, parent(x)(y)
   elseif T === U
      return true, parent(y)(x), y
   else
      return false, x, y
   end
end

function Base.promote(x::S, y::T) where {S <: RingElem, T <: RingElem}
  fl, u, v = try_promote(x, y)
  if fl
    return u, v
  else
    error("Cannot promote to common type")
  end
end

+(x::RingElem, y::RingElem) = +(promote(x, y)...)

+(x::RingElem, y::RingElement) = x + parent(x)(y)

+(x::RingElement, y::RingElem) = parent(y)(x) + y

-(x::RingElem, y::RingElem) = -(promote(x, y)...)

-(x::RingElem, y::RingElement) = x - parent(x)(y)

-(x::RingElement, y::RingElem) = parent(y)(x) - y

*(x::RingElem, y::RingElem) = *(promote(x, y)...)

*(x::RingElem, y::RingElement) = x*parent(x)(y)

*(x::RingElement, y::RingElem) = parent(y)(x)*y

function ==(x::RingElem, y::RingElem)
  fl, u, v = try_promote(x, y)
  if fl
    return u == v
  else
    return false
  end
end

==(x::RingElem, y::RingElement) = x == parent(x)(y)

==(x::RingElement, y::RingElem) = parent(y)(x) == y

include("julia/Integer.jl")

include("julia/Rational.jl")



