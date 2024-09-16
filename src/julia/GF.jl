###############################################################################
#
#   GF.jl : Julia finite fields
#
###############################################################################

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent_type(::Type{GFElem{T}}) where T <: Integer = GFField{T}

elem_type(::Type{GFField{T}}) where T <: Integer = GFElem{T}

base_ring_type(::Type{<:GFField}) = typeof(Union{})

base_ring(a::GFField) = Union{}

parent(a::GFElem) = a.parent

is_exact_type(::Type{GFElem{T}}) where T <: Integer = true

is_domain_type(::Type{GFElem{T}}) where T <: Integer = true

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::GFElem, h::UInt)
   b = 0xe08f2b4ea1cd2a13%UInt
   return xor(xor(hash(a.d), h), b)
end

@doc raw"""
    data(R::GFElem)

Return the internal data used to represent the finite field element. This
coincides with `lift` except where the internal data ids a machine integer.
"""
data(a::GFElem) = a.d

@doc raw"""
    lift(R::GFElem)

Lift the finite field element to the integers. The result will be a
multiprecision integer regardless of how the field element is represented
internally.
"""
lift(a::GFElem) = BigInt(data(a))

function zero(R::GFField{T}) where T <: Integer
   return GFElem{T}(T(0), R)
end

function one(R::GFField{T}) where T <: Integer
      return GFElem{T}(T(1), R)
end

@doc raw"""
    gen(R::GFField{T}) where T <: Integer

Return a generator of the field. Currently this returns 1.
"""
function gen(R::GFField{T}) where T <: Integer
      return GFElem{T}(T(1), R)
end

iszero(a::GFElem{T}) where T <: Integer = a.d == 0

isone(a::GFElem{T}) where T <: Integer = a.d == 1

is_unit(a::GFElem) = a.d != 0

@doc raw"""
    characteristic(R::GFField)

Return the characteristic of the given finite field.
"""
function characteristic(R::GFField)
   return R.p
end

@doc raw"""
    order(R::GFField)

Return the order, i.e. the number of element in the given finite field.
"""
function order(R::GFField)
   return R.p
end

@doc raw"""
    degree(R::GFField)

Return the degree of the given finite field.
"""
function degree(R::GFField)
   return 1
end

function deepcopy_internal(a::GFElem{T}, dict::IdDict) where T <: Integer
   R = parent(a)
   return GFElem{T}(deepcopy_internal(a.d, dict), R)
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

function expressify(x::GFElem; context = nothing)
   return expressify(x.d, context = context)
end

function show(io::IO, x::GFElem)
   print(io, x.d)
end

function show(io::IO, R::GFField)
   @show_name(io, R)
   @show_special(io, R)
   print(io, "Finite field F_", R.p)
end

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
   pmy = p - y.d
   z = x.d < pmy ? x.d + y.d : x.d - pmy
   return GFElem{T}(z, R)
end

function -(x::GFElem{T}, y::GFElem{T}) where T <: Integer
   check_parent(x, y)
   R = parent(x)
   p = characteristic(R)::T
   z = x.d < y.d ? x.d + (p - y.d) : x.d - y.d
   return GFElem{T}(z, R)
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
      y == 0 && return one(R)
      y < 0 && throw(DivideError())
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

function Base.inv(x::GFElem{T}) where T <: Integer
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

function divexact(x::GFElem{T}, y::GFElem{T}; check::Bool=true) where T <: Integer
   check_parent(x, y)
   return x*inv(y)
end

divides(a::GFElem{T}, b::GFElem{T}) where T <: Integer = true, divexact(a, b)

###############################################################################
#
#   Square root
#
###############################################################################

# helper function, computes Legendre symbol
# assumes p is an odd prime and a is in range [0, p)
legendre(a::T, p::T) where T <: Integer = powermod(a, div(p - 1, 2), p)

# Tonelli-Shanks algorithm for square root modulo a prime p
function sqrt_tonelli_shanks(a::GFElem{T}; check::Bool=true) where T <: Integer
   R = parent(a)
   p = R.p
   if p == 2 || iszero(a)
      return true, a
   end
   n = a.d
   if check && legendre(n, p) != 1
      return false, zero(R)
   end
   q = p - 1
   s = T(trailing_zeros(q))
   q = div(q, T(1) << s)
   if isone(s)
      r = powermod(n, div(p + 1, 4), p)
      return true, r
   end
   # find a quadratic nonresidue mod p
   z = T(2)
   while z < p
      if legendre(z, p) == p - 1
         break
      end
      z += 1
   end
   c = powermod(z, q, p)
   r = powermod(n, div(q + 1, 2), p)
   t = powermod(n, q, p)
   m = s
   t2 = zero(R)
   while !isone(t)
      t2 = mulmod(t, t, p)
      i = T(1)
      while i <= m
         if isone(t2)
            break
         end
         t2 = mulmod(t2, t2, p)
         i += 1
      end
      b = powermod(c, 1 << (m - i - 1), p)
      r = mulmod(r, b, p)
      c = mulmod(b, b, p)
      t = mulmod(t, c, p)
      m = i
   end
   return true, r
end

function Base.sqrt(a::GFElem{T}; check::Bool=true) where T <: Integer
   f1, s1 = sqrt_tonelli_shanks(a; check=check)
   check && !f1 && error("Not a square in sqrt")
   return s1
end

function is_square(a::GFElem{T}) where T <: Integer
   f1, s1 = sqrt_tonelli_shanks(a; check=true)
   return f1
end

function is_square_with_sqrt(a::GFElem{T}) where T <: Integer
   return sqrt_tonelli_shanks(a; check=true)
end

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

function add!(z::GFElem{BigInt}, x::GFElem{BigInt})
   R = parent(x)
   p = R.p::BigInt
   d = add!(z.d, x.d)
   if d < p
      return GFElem{BigInt}(d, R)
   else
      return GFElem{BigInt}(d - p, R)
   end
end

function add!(z::GFElem{BigInt}, x::GFElem{BigInt}, y::GFElem{BigInt})
   R = parent(x)
   p = R.p::BigInt
   d = add!(z.d, x.d, y.d)
   if d < p
      return GFElem{BigInt}(d, R)
   else
      return GFElem{BigInt}(d - p, R)
   end
end

###############################################################################
#
#   Random functions
#
###############################################################################

Random.Sampler(RNG::Type{<:AbstractRNG}, R::FinField, n::Random.Repetition) =
  Random.SamplerSimple(R, Random.Sampler(RNG, BigInt(0):BigInt(characteristic(R) - 1), n))

rand(rng::AbstractRNG, R::Random.SamplerSimple{GFField{T}}) where T =
   GFElem{T}(rand(rng, R.data), R[])

function rand(rng::AbstractRNG, Ksp::Random.SamplerSimple{<:FinField})
   K = Ksp[]
   r = degree(K)
   alpha = gen(K)
   res = zero(K)
   for i = 0:(r-1)
      c = rand(rng, Ksp.data)
      res += c * alpha^i
   end
   return res
end

Random.gentype(T::Type{<:FinField}) = elem_type(T)

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

@doc raw"""
    GF(p::T; check::Bool=true) where T <: Integer

Return the finite field $\mathbb{F}_p$, where $p$ is a prime.
By default, the integer $p$ is checked with a probabilistic algorithm for primality.
When `check == false`, no check is made, but the behaviour of the resulting object
is undefined if $p$ is composite.
"""
function GF(p::T; cached::Bool = true, check::Bool=true) where T <: Integer
   check && !is_probable_prime(p) && throw(DomainError(p, "Characteristic is not prime in GF(p)"))
   return GFField{T}(p; cached = cached)
end
