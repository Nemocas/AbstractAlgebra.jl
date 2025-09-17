################################################################################
#
#   Fields.jl : generic fields
#
################################################################################

is_domain_type(::Type{T}) where {T <: FieldElem} = true

is_zero_divisor(a::T) where T <: FieldElem = is_zero(a)

is_unit(a::FieldElem) = !iszero(a)

//(a::T, b::T) where {T <: FieldElem} = divexact(a, b)

//(x::T, y::Union{Integer, Rational}) where {T <: RingElem} = x//parent(x)(y)
                                          
//(x::Union{Integer, Rational}, y::T) where {T <: RingElem} = parent(y)(x)//y

Base.divrem(a::T, b::T) where {T <: FieldElem} = divexact(a, b), zero(parent(a))

Base.div(a::T, b::T) where {T <: FieldElem} = divexact(a, b)

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

function factor(x::FieldElement)
  @req !is_zero(x) "Element must be non-zero"
  return Fac(x, Dict{typeof(x), Int}())
end

function canonical_unit(x::FieldElement)
  iszero(x) && return one(x)
  return x
end

@doc raw"""
   exercise_function(n::Int)

Return n+1 given the number n
A Practice for the Summer School Exercise

# Examples:
```jldoctest 
julia> exercise_function(20)
21

julia> AbstractAlgebra.exercise_function(10)
11
```
"""
function exercise_function(n::Int)
   # create an integer 10^(ceil(n/3))
   #  k = ceil(Int, n / 3 + (n % 3 != 0))
   #  b_int = big(10)^k
   #  # turn into a rational in rational numbers 
   #  b = QQ(b_int)

   #  # define polynomial ring over rational numbers
   #  Qx, x = QQ["x"]

   #  # define the cubic polynomial f = x^3 + 3*b*x + 3
   #  f = x^3 + 3*b*x + 3

   #  #define the number field
   #  #K, _ = number_field(f)

    return n+1
end
