################################################################################
#
#   Fields.jl : generic fields
#
################################################################################

include("julia/GF.jl")

//(a::T, b::T) where {T <: FieldElem} = divexact(a, b)

//(x::T, y::Union{Integer, Rational}) where {T <: RingElem} = x//parent(x)(y)
                                          
//(x::Union{Integer, Rational}, y::T) where {T <: RingElem} = parent(y)(x)//y

Base.divrem(a::T, b::T) where {T <: FieldElem} = divexact(a, b), zero(parent(a))

div(a::T, b::T) where {T <: FieldElem} = divexact(a, b)

function gcd(x::T, y::T) where {T <: FieldElem}
   check_parent(x, y)
   return iszero(x) && iszero(y) ? zero(parent(y)) : one(parent(y))
end

function gcdx(x::T, y::T) where {T <: FieldElem}
   check_parent(x, y)
   R = parent(x)
   if iszero(x)
      if iszero(y)
         return zero(R), one(R), zero(R)
      end
      return one(R), zero(R), inv(y)
   else
      return one(R), inv(x), zero(R)
   end
end

characteristic(R::Field) = 0


