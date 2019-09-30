###############################################################################
#
#   GF.jl : Julia finite fields
#
###############################################################################

export GF

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent_type(::Type{GFElem{T}}) where T <: Integer = GFField{T}

elem_type(::Type{GFField{T}}) where T <: Integer = GFElem{T}

@doc Markdown.doc"""
    base_ring(a::GFField)
> Return `Union{}` as this field is not dependent on another field.
"""
base_ring(a::GFField) = Union{}

@doc Markdown.doc"""
    base_ring(a::GFElem)
> Return `Union{}` as this field is not dependent on another field.
"""
base_ring(a::GFElem) = Union{}

@doc Markdown.doc"""
    parent(a::GFElem)
> Return the parent of the given finite field element.
"""
parent(a::GFElem) = a.parent

isexact_type(::Type{GFElem{T}}) where T <: Integer = true

isdomain_type(::Type{GFElem{T}}) where T <: Integer = true

function check_parent(a::GFElem, b::GFElem)
   a.parent != b.parent && error("Operations on distinct finite fields not supported")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::GFElem, h::UInt)
   b = 0xe08f2b4ea1cd2a13%UInt
   return xor(xor(hash(a.d), h), b)
end

@doc Markdown.doc"""
    zero(R::GFField{T}) where T <: Integer
> Return the additive identity, zero, in the given finite field.
"""
function zero(R::GFField{T}) where T <: Integer
   return GFElem{T}(T(0), R)
end

@doc Markdown.doc"""
    one(R::GFField{T}) where T <: Integer
> Return the additive identity, zero, in the given finite field.
"""
function one(R::GFField{T}) where T <: Integer
      return GFElem{T}(T(1), R)
end

@doc Markdown.doc"""
    gen(R::GFField{T}) where T <: Integer
> Return a generator of the field. Currently this returns 1.
"""
function gen(R::GFField{T}) where T <: Integer
      return GFElem{T}(T(1), R)
end

@doc Markdown.doc"""
    iszero(a::GFElem{T}) where T <: Integer
> Return true if the given element of the finite field is zero.
"""
iszero(a::GFElem{T}) where T <: Integer = a.d == 0

@doc Markdown.doc"""
    isone(a::GFElem{T}) where T <: Integer
> Return true if the given element of the finite field is one.
"""
isone(a::GFElem{T}) where T <: Integer = a.d == 1

@doc Markdown.doc"""
    isunit(a::GFElem)
> Return `true` if the given finite field element is invertible, i.e. nonzero,
> otherwise return `false`.
"""
isunit(a::GFElem) = a.d != 0

@doc Markdown.doc"""
    characteristic(R::GFField)
> Return the characteristic of the given finite field.
"""
function characteristic(R::GFField)
   return R.p
end

@doc Markdown.doc"""
    order(R::GFField)
> Return the order, i.e. the number of element in the given finite field.
"""
function order(R::GFField)
   return R.p
end

@doc Markdown.doc"""
    degree(R::GFField)
> Return the degree of the given finite field.
"""
function degree(R::GFField)
   return 1
end

function deepcopy_internal(a::GFElem{T}, dict::IdDict) where T <: Integer
   R = parent(a)
   return GFElem{T}(deepcopy(a.d), R)
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::GFElem) = x

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, x::GFElem)
   print(io, x.d)
end

function show(io::IO, R::GFField)
   print(io, "Finite field F_", R.p)
end

needs_parentheses(x::GFElem) = false

displayed_with_minus_in_front(x::GFElem) = false

show_minus_one(::Type{GFElem{T}}) where T <: Integer = true

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(x::GFElem{T}) where T <: Integer
   if x.d == 0
      return deepcopy(x)
   else
      R = parent(x)
      return GFElem{T}(R.p - x.d, R)
   end
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(x::GFElem{T}, y::GFElem{T}) where T <: Integer
   check_parent(x, y)
   R = parent(x)
   p = characteristic(R)::T
   d = x.d + y.d - p
   if d < 0
      return GFElem{T}(d + p, R)
   else
      return GFElem{T}(d, R)
   end
end

function -(x::GFElem{T}, y::GFElem{T}) where T <: Integer
   check_parent(x, y)
   R = parent(x)
   p = characteristic(R)::T
   d = x.d - y.d
   if d < 0
      return GFElem{T}(d + p, R)
   else
      return GFElem{T}(d, R)
   end
end

function *(x::GFElem{T}, y::GFElem{T}) where T <: Integer
   check_parent(x, y)
   R = parent(x)
   return R(widen(x.d)*widen(y.d))
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::Integer, y::GFElem{T}) where T <: Integer
   R = parent(y)
   return R(widen(x)*widen(y.d))
end

*(x::GFElem{T}, y::Integer) where T <: Integer = y*x

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::GFElem{T}, y::Integer) where T <: Integer
   R = parent(x)
   p = R.p::T
   if x.d == 0
      return deepcopy(x)
   end
   if y < 0
      x = inv(x)
      y = -y
   end
   if y >= p - 1
      y1 = y%(p - 1)
      y = convert(T, y1)::T
   else
      y = convert(T, y)::T
   end
   if y == 0
      return one(R)
   elseif y == 1
      return deepcopy(x)
   end
   bit = T(1) << (ndigits(y, base = 2) - 1)
   z = x
   bit >>= 1
   while bit != 0
      z = z*z
      if (bit & y) != 0
         z *= x
      end
      bit >>= 1
   end
   return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::GFElem{T}, y::GFElem{T}) where T <: Integer
   check_parent(x, y)
   return x.d == y.d
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(x::GFElem{T}) where T <: Integer
   x == 0 && throw(DivideError())
   R = parent(x)
   p = R.p::T
   g, s, t = gcdx(x.d, p)
   g != 1 && error("Characteristic not prime in ", R)
   return R(s)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::GFElem{T}, y::GFElem{T}) where T <: Integer
   check_parent(x, y)
   return x*inv(y)
end

divides(a::GFElem{T}, b::GFElem{T}) where T <: Integer = true, divexact(a, b)

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(z::GFElem{T}) where T <: Integer
   R = parent(z)
   d = zero!(z.d)
   return GFElem{T}(d, R)
end

function mul!(z::GFElem{T}, x::GFElem{T}, y::GFElem{T}) where T <: Integer
   return x*y
end

function mul!(z::GFElem{BigInt}, x::GFElem{BigInt}, y::GFElem{BigInt})
   R = parent(x)
   p = R.p::BigInt
   d = mul!(z.d, x.d, y.d)
   if d >= p
      return GFElem{BigInt}(d%p, R)
   else
      return GFElem{BigInt}(d, R)
   end
end

function addeq!(z::GFElem{T}, x::GFElem{T}) where T <: Integer
   R = parent(x)
   p = R.p::T
   d = addeq!(z.d, x.d)
   if d < p
      return GFElem{T}(d, R)
   else
      return GFElem{T}(d - p, R)
   end
end

function add!(z::GFElem{T}, x::GFElem{T}, y::GFElem{T}) where T <: Integer
   R = parent(x)
   p = R.p::T
   d = add!(z.d, x.d, y.d)
   if d < p
      return GFElem{T}(d, R)
   else
      return GFElem{T}(d - p, R)
   end
end

###############################################################################
#
#   Random functions
#
###############################################################################

Random.Sampler(RNG::Type{<:AbstractRNG}, R::GFField, n::Random.Repetition) =
   Random.SamplerSimple(R, Random.Sampler(RNG, 0:R.p - 1, n))

rand(rng::AbstractRNG, R::Random.SamplerSimple{GFField{T}}) where T =
   GFElem{T}(rand(rng, R.data), R[])

Random.gentype(T::Type{<:GFField}) = elem_type(T)

###############################################################################
#
#   Promotions
#
###############################################################################

promote_rule(::Type{GFElem{S}}, ::Type{T}) where {S <: Integer, T <: Integer} = GFElem{S}

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::GFField{T})() where T <: Integer
   return GFElem{T}(T(0), R)
end

function (R::GFField{T})(a::Integer) where T <: Integer
   p = R.p::T
   d = convert(T, a%p)::T
   if d < 0
      d += p
   end
   return GFElem{T}(d, R)
end

function (R::GFField{T})(a::GFElem{T}) where T <: Integer
   parent(a) != R && error("Coercion between finite fields not implemented")
   return a
end

###############################################################################
#
#   GF(p) constructor
#
###############################################################################

@doc Markdown.doc"""
    GF(p::T) where T <: Integer
> Return the finite field $\mathbb{F}_p$, where $p$ is a prime. The integer
> $p$ is not checked for primality, but the behaviour of the resulting object
> is undefined if $p$ is composite.
"""
function GF(p::T) where T <: Integer
   p <= 0 && error("Characteristic is not prime in GF(p)")
   return GFField{T}(p)
end
