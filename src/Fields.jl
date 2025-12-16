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
    evaluation_points(K::Field, n::Int) -> Vector{FieldElem}

If $K$ contains at least $n$ elements, this function returns a vector $v$ with $n$ distinct elements from $K$. Otherwise it returns an empty vector. 
"""
function evaluation_points(K::Field, n::Int)
   @req n > 0 "n must be positive."
   if is_finite(K)
      v = elem_type(K)[]
      if order(K) < n
         return v
      else
         append!(v, Iterators.take(K, n))
      end
   else
      @assert is_zero(characteristic(K)) "This should not happen. Please open a bug report."
      v = Vector{elem_type(K)}(undef, n)
      v[1] = K(div(-n, 2))
      for i = 2:n
         v[i] = v[i-1] + 1
      end
   end
   return v
end

function evaluation_points(K::Generic.RationalFunctionField, n::Int)
   @req n > 0 "n must be positive."
   F = base_ring(K)
   v = evaluation_points(F, n)
   if is_empty(v)
      v = elem_type(K)[]
      base_v = evaluation_points(F, Int(order(F)))
      d = clog(order(F), ZZ(n))
      genKpowers = powers(gen(K), d - 1)
      for coeffs in ProductIterator(base_v, d; inplace=true)
         push!(v, sum(coeffs[j] * genKpowers[j] for j in 1:d))
         if length(v) == n
            break
         end
      end
   end
   return v
end

function ext_of_degree(K::Field, n::Int)
   error("Cannot extend field")
end
