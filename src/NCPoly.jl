###############################################################################
#
#   NCPoly.jl : polynomials over noncommutative rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{NCPolyRing{T}}) where T <: NCRingElem = parent_type(T)

base_ring(R::NCPolyRing{T}) where T <: NCRingElem = R.base_ring::parent_type(T)

coefficient_ring(R::NCPolyRing) = base_ring(R)

function is_exact_type(a::Type{T}) where {S <: NCRingElem, T <: NCPolyRingElem{S}}
   return is_exact_type(S)
end

@doc raw"""
    dense_poly_type(::Type{T}) where T<:NCRingElement
    dense_poly_type(::T) where T<:NCRingElement
    dense_poly_type(::Type{S}) where S<:NCRing
    dense_poly_type(::S) where S<:NCRing

The type of univariate polynomials with coefficients of type `T` respectively `elem_type(S)`.
Falls back to `Generic.NCPoly{T}` respectively `Generic.Poly{T}`.

See also [`dense_poly_ring_type`](@ref), [`mpoly_type`](@ref) and [`mpoly_ring_type`](@ref).

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> dense_poly_type(AbstractAlgebra.ZZ(1))
AbstractAlgebra.Generic.Poly{BigInt}

julia> dense_poly_type(elem_type(AbstractAlgebra.ZZ))
AbstractAlgebra.Generic.Poly{BigInt}

julia> dense_poly_type(AbstractAlgebra.ZZ)
AbstractAlgebra.Generic.Poly{BigInt}

julia> dense_poly_type(typeof(AbstractAlgebra.ZZ))
AbstractAlgebra.Generic.Poly{BigInt}
```
"""
dense_poly_type(::Type{T}) where T<:NCRingElement = Generic.NCPoly{T}
dense_poly_type(::Type{S}) where S<:NCRing = dense_poly_type(elem_type(S))
dense_poly_type(x) = dense_poly_type(typeof(x)) # to stop this method from eternally recursing on itself, we better add ...
dense_poly_type(::Type{T}) where T = throw(ArgumentError("Type `$T` must be subtype of `NCRingElement`."))

@doc raw"""
    dense_poly_ring_type(::Type{T}) where T<:NCRingElement
    dense_poly_ring_type(::T) where T<:NCRingElement
    dense_poly_ring_type(::Type{S}) where S<:NCRing
    dense_poly_ring_type(::S) where S<:NCRing

The type of univariate polynomial rings with coefficients of type `T` respectively
`elem_type(S)`. Implemented via [`dense_poly_type`](@ref).

See also [`mpoly_type`](@ref) and [`mpoly_ring_type`](@ref).

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> dense_poly_ring_type(AbstractAlgebra.ZZ(1))
AbstractAlgebra.Generic.PolyRing{BigInt}

julia> dense_poly_ring_type(elem_type(AbstractAlgebra.ZZ))
AbstractAlgebra.Generic.PolyRing{BigInt}

julia> dense_poly_ring_type(AbstractAlgebra.ZZ)
AbstractAlgebra.Generic.PolyRing{BigInt}

julia> dense_poly_ring_type(typeof(AbstractAlgebra.ZZ))
AbstractAlgebra.Generic.PolyRing{BigInt}
```
"""
dense_poly_ring_type(x) = parent_type(dense_poly_type(x))

@doc raw"""
    var(a::NCPolyRing)

Return the internal name of the generator of the polynomial ring. Note that
this is returned as a `Symbol` not a `String`.
"""
var(a::NCPolyRing) = a.S

@doc raw"""
    symbols(a::NCPolyRing)

Return an array of the variable names for the polynomial ring. Note that
this is returned as an array of `Symbol` not `String`.
"""
symbols(a::NCPolyRing) = [a.S]

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::NCPolyRingElem, h::UInt)
   b = 0xd3f41ffbf953cbd8%UInt
   for i in 0:length(a) - 1
      b = xor(b, xor(hash(coeff(a, i), h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

zero(R::NCPolyRing) = R(0)

one(R::NCPolyRing) = R(1)

@doc raw"""
    gen(R::NCPolyRing)

Return the generator of the given polynomial ring.
"""
gen(R::NCPolyRing) = R([zero(base_ring(R)), one(base_ring(R))])

is_term(a::T) where T <: NCRingElem = true

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, p::NCPolyRing)
   io = pretty(io)
   print(io, "Univariate polynomial ring in ", var(p), " over ")
   print(IOContext(io, :compact => true), Lowercase(), base_ring(p))
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(a::NCPolyRingElem{T}, b::NCPolyRingElem{T}) where T <: NCRingElem
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) + coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, deepcopy(coeff(b, i - 1)))
      i += 1
   end
   z = set_length!(z, normalise(z, i - 1))
   return z
end

function -(a::NCPolyRingElem{T}, b::NCPolyRingElem{T}) where T <: NCRingElem
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) - coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, -coeff(b, i - 1))
      i += 1
   end
   z = set_length!(z, normalise(z, i - 1))
   return z
end

function *(a::NCPolyRingElem{T}, b::NCPolyRingElem{T}) where T <: NCRingElem
   lena = length(a)
   lenb = length(b)
   if lena == 0 || lenb == 0
      return parent(a)()
   end
   t = base_ring(a)()
   lenz = lena + lenb - 1
   d = Vector{T}(undef, lenz)
   for i = 1:lena
      d[i] = coeff(a, i - 1)*coeff(b, 0)
   end
   for i = 2:lenb
      d[lena + i - 1] = a.coeffs[lena]*coeff(b, i - 1)
   end
   for i = 1:lena - 1
      for j = 2:lenb
         t = mul!(t, coeff(a, i - 1), b.coeffs[j])
         d[i + j - 1] = addeq!(d[i + j - 1], t)
      end
   end
   z = parent(a)(d)
   z = set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::NCPolyRingElem{T}) where T <: NCRingElem
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   return z
end

function *(a::NCPolyRingElem{T}, b::T) where T <: NCRingElem
   len = length(a)
   z = parent(a)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, coeff(a, i - 1)*b)
   end
   z = set_length!(z, normalise(z, len))
   return z
end

function +(a::T, b::NCPolyRingElem{T}) where {T <: NCRingElem}
   z = deepcopy(b)
   len = length(z)
   z = setcoeff!(z, 0, a + coeff(b, 0))
   z = set_length!(z, normalise(z, len))
   return z
end

+(a::NCPolyRingElem{T}, b::T) where {T <: NCRingElem} =  b + a

+(a::Union{Integer, Rational}, b::NCPolyRingElem{T}) where {T <: NCRingElem} =  parent(b)(a) + b

+(a::NCPolyRingElem{T}, b::Union{Integer, Rational}) where {T <: NCRingElem} =  b + a

-(a::Union{Integer, Rational}, b::NCPolyRingElem{T}) where {T <: NCRingElem} =  parent(b)(a) - b

-(a::NCPolyRingElem{T}, b::Union{Integer, Rational}) where {T <: NCRingElem} =  a - parent(a)(b)

###############################################################################
#
#   Powering
#
###############################################################################

@doc raw"""
    ^(a::NCPolyRingElem{T}, b::Int) where T <: NCRingElem

Return $a^b$. We require $b \geq 0$.
"""
function ^(a::NCPolyRingElem{T}, b::Int) where T <: NCRingElem
   b < 0 && throw(DomainError(b, "exponent must be >= 0"))
   # special case powers of x for constructing polynomials efficiently
   R = parent(a)
   if is_gen(a)
      z = R()
      fit!(z, b + 1)
      z = setcoeff!(z, b, deepcopy(coeff(a, 1)))
      for i = 1:b
         z = setcoeff!(z, i - 1, deepcopy(coeff(a, 0)))
      end
      z = set_length!(z, b + 1)
      return z
   elseif b == 0
      return one(R)
   elseif length(a) == 0
      return zero(R)
   elseif length(a) == 1
      return R(coeff(a, 0)^b)
   elseif b == 1
      return deepcopy(a)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit != 0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

###############################################################################
#
#   Comparisons
#
###############################################################################

@doc raw"""
    ==(x::NCPolyRingElem{T}, y::NCPolyRingElem{T}) where T <: NCRingElem

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(x::NCPolyRingElem{T}, y::NCPolyRingElem{T}) where T <: NCRingElem
   b = check_parent(x, y, false)
   !b && return false
   if length(x) != length(y)
      return false
   else
      for i = 1:length(x)
         if coeff(x, i - 1) != coeff(y, i - 1)
            return false
         end
      end
   end
   return true
end

@doc raw"""
    isequal(x::NCPolyRingElem{T}, y::NCPolyRingElem{T}) where T <: NCRingElem

Return `true` if $x == y$ exactly, otherwise return `false`. This function is
useful in cases where the coefficients of the polynomial are inexact, e.g.
power series. Only if the power series are precisely the same, to the same
precision, are they declared equal by this function.
"""
function isequal(x::NCPolyRingElem{T}, y::NCPolyRingElem{T}) where T <: NCRingElem
   if parent(x) != parent(y)
      return false
   end
   if length(x) != length(y)
      return false
   end
   for i = 1:length(x)
      if !isequal(coeff(x, i - 1), coeff(y, i - 1))
         return false
      end
   end
   return true
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

@doc raw"""
    ==(x::NCPolyRingElem{T}, y::T) where T <: NCRingElem

Return `true` if $x == y$.
"""
==(x::NCPolyRingElem{T}, y::T) where T <: NCRingElem = ((length(x) == 0 && y == 0)
                        || (length(x) == 1 && coeff(x, 0) == y))

@doc raw"""
    ==(x::T, y::NCPolyRingElem{T}) where T <: NCRingElem

Return `true` if $x = y$.
"""
==(x::T, y::NCPolyRingElem{T}) where T <: NCRingElem = y == x

@doc raw"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::NCPolyRingElem)

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::NCPolyRingElem) = y == x

###############################################################################
#
#   Truncation
#
###############################################################################

@doc raw"""
    mullow(a::NCPolyRingElem{T}, b::NCPolyRingElem{T}, n::Int) where T <: NCRingElem

Return $a\times b$ truncated to $n$ terms.
"""
function mullow(a::NCPolyRingElem{T}, b::NCPolyRingElem{T}, n::Int) where T <: NCRingElem
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   if lena == 0 || lenb == 0
      return zero(parent(a))
   end
   if n < 0
      n = 0
   end
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, n)
   d = Vector{T}(undef, lenz)
   for i = 1:min(lena, lenz)
      d[i] = coeff(a, i - 1)*coeff(b, 0)
   end
   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = coeff(a, lena - 1)*coeff(b, j - 1)
      end
   end
   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            t = mul!(t, coeff(a, i - 1), coeff(b, j - 1))
            d[i + j - 1] = addeq!(d[i + j - 1], t)
         end
      end
   end
   z = parent(a)(d)
   z = set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

@doc raw"""
    divexact_right(f::NCPolyRingElem{T}, g::NCPolyRingElem{T}; check::Bool=true) where T <: NCRingElem

Assuming $f = qg$, return $q$. By default if the division is not exact an
exception is raised. If `check=false` this test is omitted.
"""
function divexact_right(f::NCPolyRingElem{T}, g::NCPolyRingElem{T}; check::Bool=true) where T <: NCRingElem
   check_parent(f, g)
   iszero(g) && throw(DivideError())
   if iszero(f)
      return zero(parent(f))
   end
   lenq = length(f) - length(g) + 1
   d = Vector{T}(undef, lenq)
   for i = 1:lenq
      d[i] = zero(base_ring(f))
   end
   x = gen(parent(f))
   leng = length(g)
   while length(f) >= leng
      lenf = length(f)
      q1 = d[lenf - leng + 1] = divexact_right(coeff(f, lenf - 1), coeff(g, leng - 1); check=check)
      f = f - shift_left(q1*g, lenf - leng)
      if length(f) == lenf # inexact case
         f = set_length!(f, normalise(f, lenf - 1))
      end
   end
   check && length(f) != 0 && error("Not an exact division")
   q = parent(f)(d)
   q = set_length!(q, lenq)
   return q
end

@doc raw"""
    divexact_left(f::NCPolyRingElem{T}, g::NCPolyRingElem{T}; check::Bool=true) where T <: NCRingElem

Assuming $f = gq$, return $q$. By default if the division is not exact an
exception is raised. If `check=false` this test is omitted.
"""
function divexact_left(f::NCPolyRingElem{T}, g::NCPolyRingElem{T}; check::Bool=true) where T <: NCRingElem
   check_parent(f, g)
   iszero(g) && throw(DivideError())
   if iszero(f)
      return zero(parent(f))
   end
   lenq = length(f) - length(g) + 1
   d = Vector{T}(undef, lenq)
   for i = 1:lenq
      d[i] = zero(base_ring(f))
   end
   x = gen(parent(f))
   leng = length(g)
   while length(f) >= leng
      lenf = length(f)
      q1 = d[lenf - leng + 1] = divexact_left(coeff(f, lenf - 1), coeff(g, leng - 1); check=check)
      f = f - shift_left(g*q1, lenf - leng)
      if length(f) == lenf # inexact case
         f = set_length!(f, normalise(f, lenf - 1))
      end
   end
   check && length(f) != 0 && error("Not an exact division")
   q = parent(f)(d)
   q = set_length!(q, lenq)
   return q
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

@doc raw"""
    divexact_right(a::NCPolyRingElem{T}, b::T; check::Bool=true) where T <: NCRingElem

Assuming $a = qb$, return $q$. By default if the division is not exact an
exception is raised. If `check=false` this test is omitted.
"""
function divexact_right(a::NCPolyRingElem{T}, b::T; check::Bool=true) where T <: NCRingElem
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact_right(coeff(a, i - 1), b; check=check))
   end
   z = set_length!(z, length(a))
   return z
end

@doc raw"""
    divexact_left(a::NCPolyRingElem{T}, b::T; check::Bool=true) where T <: NCRingElem

Assuming $a = bq$, return $q$. By default if the division is not exact an
exception is raised. If `check=false` this test is omitted.
"""
function divexact_left(a::NCPolyRingElem{T}, b::T; check::Bool=true) where T <: NCRingElem
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact_left(coeff(a, i - 1), b; check=check))
   end
   z = set_length!(z, length(a))
   return z
end

@doc raw"""
    divexact_right(a::NCPolyRingElem, b::Union{Integer, Rational, AbstractFloat}; check::Bool=true)

Assuming $a = qb$, return $q$. By default if the division is not exact an
exception is raised. If `check=false` this test is omitted.
"""
function divexact_right(a::NCPolyRingElem, b::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact_right(coeff(a, i - 1), b; check=check))
   end
   z = set_length!(z, length(a))
   return z
end

@doc raw"""
    divexact_left(a::NCPolyRingElem, b::Union{Integer, Rational, AbstractFloat}; check::Bool=true)

Assuming $a = bq$, return $q$. By default if the division is not exact an
exception is raised. If `check=false` this test is omitted.
"""
function divexact_left(a::NCPolyRingElem, b::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   return divexact_right(a, b; check=check)
end

###############################################################################
#
#   Evaluation
#
###############################################################################

@doc raw"""
    evaluate(a::NCPolyRingElem, b::T) where T <: NCRingElem

Evaluate the polynomial $a$ at the value $b$ and return the result.
"""
function evaluate(a::NCPolyRingElem, b::T) where T <: NCRingElem
   i = length(a)
   R = base_ring(a)
   if i == 0
       return zero(R)
   end
   z = R(coeff(a, i - 1))
   while i > 1
      i -= 1
      z = R(coeff(a, i - 1)) + z*b
      parent(z) # To work around a bug in julia
   end
   return z
end

@doc raw"""
    evaluate(a::NCPolyRingElem, b::Union{Integer, Rational, AbstractFloat})

Evaluate the polynomial $a$ at the value $b$ and return the result.
"""
function evaluate(a::NCPolyRingElem, b::Union{Integer, Rational, AbstractFloat})
   i = length(a)
   R = base_ring(a)
   if i == 0
       return zero(R)
   end
   z = R(coeff(a, i - 1))
   while i > 1
      i -= 1
      z = R(coeff(a, i - 1)) + z*b
      parent(z) # To work around a bug in julia
   end
   return z
end

# Note: composition is not associative, e.g. consider fo(goh) vs (fog)oh
# for f and g of degree 2 and h of degree 1 -- and recall coeffs don't commute

################################################################################
#
#  Change base ring
#
################################################################################

function change_base_ring(R::NCRing, p::NCPolyRingElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p), cached)) where T <: NCRingElement
   return _map(R, p, parent)
end

function change_coefficient_ring(R::NCRing, p::NCPolyRingElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p), cached)) where T <: NCRingElement
  return change_base_ring(R, p; cached = cached, parent = parent)
end

################################################################################
#
#  Map
#
################################################################################

_make_parent(g::T, p::NCPolyRingElem, cached::Bool) where {T} =
   _change_poly_ring(parent(g(zero(base_ring(p)))),
                     parent(p), cached)

function map_coefficients(g::T, p::NCPolyRingElem{<:NCRingElement};
                    cached::Bool = true,
		    parent::NCPolyRing = _make_parent(g, p, cached)) where {T}
   return _map(g, p, parent)
end

function _map(g::T, p::NCPolyRingElem, Rx) where {T}
   R = base_ring(Rx)
   new_coefficients = elem_type(R)[let c = coeff(p, i)
                                     iszero(c) ? zero(R) : R(g(c))
                                   end for i in 0:degree(p)]
   return Rx(new_coefficients)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addmul!(z::NCPolyRingElem{T}, x::NCPolyRingElem{T}, y::NCPolyRingElem{T}, c::NCPolyRingElem{T}) where T <: NCRingElem
   c = mul!(c, x, y)
   z = addeq!(z, c)
   return z
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::NCPolyRing, dr::AbstractUnitRange{Int}, _) = elem_type(S)

function RandomExtensions.make(S::NCPolyRing, deg_range::AbstractUnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, deg_range, vs[1]) # forward to default Make constructor
   else
      Make(S, deg_range, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make3{<:NCPolyRingElem,
                                         <:NCPolyRing,
                                         <:AbstractUnitRange{Int}}})
   S, deg_range, v = sp[][1:end]
   R = base_ring(S)
   f = S()
   x = gen(S)
   for i = 0:rand(rng, deg_range)
      f += rand(rng, v)*x^i
   end
   return f
end

rand(rng::AbstractRNG, S::NCPolyRing, deg_range::AbstractUnitRange{Int}, v...) =
   rand(rng, make(S, deg_range, v...))

rand(S::NCPolyRing, deg_range, v...) = rand(Random.GLOBAL_RNG, S, deg_range, v...)


###############################################################################
#
#   polynomial_ring constructor
#
###############################################################################

@doc raw"""
    polynomial_ring(R::NCRing, s::VarName = :x; cached::Bool = true)

Given a base ring `R` and symbol/string `s` specifying how the generator
(variable) should be printed, return a tuple `S, x` representing the new
polynomial ring $S = R[x]$ and the generator $x$ of the ring.

By default the parent object `S` depends only on `R` and `x` and will be cached.
Setting the optional argument `cached` to `false` will prevent the parent object `S` from being cached.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over univariate polynomial ring, y)
```
"""
function polynomial_ring(R::NCRing, s::VarName =:x; kw...)
   S = polynomial_ring_only(R, Symbol(s); kw...)
   (S, gen(S))
end

@doc raw"""
    polynomial_ring_only(R::NCRing, s::Symbol; cached::Bool=true)

Like [`polynomial_ring(R::NCRing, s::Symbol)`](@ref) but return only the
polynomial ring.
"""
polynomial_ring_only(R::T, s::Symbol; cached::Bool=true) where T<:NCRing =
   dense_poly_ring_type(T)(R, s, cached)

# Simplified constructor

PolyRing(R::NCRing) = polynomial_ring_only(R, :x; cached=false)
