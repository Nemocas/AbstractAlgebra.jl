###############################################################################
#
#   Residue.jl : residue rings (modulo a principal ideal)
#
###############################################################################

export ResidueRing, add!, addeq!, data, modulus, mul!, lift 

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring(S::ResRing{T}) where {T <: RingElement} = S.base_ring::parent_type(T)

base_ring(r::ResElem) = base_ring(parent(r))

parent(a::ResElem) = a.parent

isdomain_type(a::Type{T}) where T <: ResElem = false

function isexact_type(a::Type{T}) where {S <: RingElement, T <: ResElem{S}}
   return isexact_type(S)
end

function check_parent_type(a::ResRing{T}, b::ResRing{T}) where {T <: RingElement}
   # exists only to check types of parents agree
end

function check_parent(a::ResElem, b::ResElem, throw::Bool = true)
   if parent(a) != parent(b)
      check_parent_type(parent(a), parent(b))
      fl = modulus(parent(a)) != modulus(parent(b))
      fl && throw && error("Incompatible moduli in residue operation")
      return !fl
      #CF: maybe extend to divisibility?
   end
   return true
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::ResElem, h::UInt)
   b = 0x539c1c8715c1adc2%UInt
   return xor(b, xor(hash(data(a), h), h))
end

@doc Markdown.doc"""
    modulus(R::ResRing)

Return the modulus $a$ of the given residue ring $S = R/(a)$.
"""
function modulus(S::ResRing)
   return S.modulus
end

@doc Markdown.doc"""
    modulus(R::ResElem)

Return the modulus $a$ of the residue ring $S = R/(a)$ that the supplied
residue $r$ belongs to.
"""
function modulus(r::ResElem)
   return modulus(parent(r))
end

data(a::ResElem) = a.data

lift(a::ResElem) = data(a)

lift(a::ResElem{Int}) = BigInt(data(a))

zero(R::ResRing) = R(0)

one(R::ResRing) = R(1)

iszero(a::ResElem) = iszero(data(a))

isone(a::ResElem) = isone(data(a))

function isunit(a::ResElem)
   g = gcd(data(a), modulus(a))
   return isone(g)
end

deepcopy_internal(a::ResElem, dict::IdDict) =
   parent(a)(deepcopy(data(a)))

function characteristic(a::ResRing{T}) where T <: Integer
   return modulus(a)
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

function canonical_unit(x::ResElem{<:Union{Integer, RingElem}})
 #the simple return x does not work
  # - if x == 0, this is not a unit
  # - if R is not a field....
  if iszero(x)
    return one(parent(x))
  end
  g = gcd(modulus(x), data(x))
  u = divexact(data(x), g)
  a, b = ppio(modulus(x), u)
  if isone(a)
    r = u
  elseif isone(b)
    r = b
  else
    r = crt(one(parent(a)), a, u, b)
  end
  return parent(x)(r)
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(a::ResElem; context = nothing)
   return expressify(data(a), context = context)
end

function show(io::IO, ::MIME"text/plain", a::ResElem)
  print(io, obj_to_string(a, context = io))
end

function show(io::IO, a::ResElem)
  print(io, obj_to_string(a, context = io))
end

function show(io::IO, a::ResRing)
   print(IOContext(io, :compact => true), "Residue ring of ", base_ring(a), " modulo ", modulus(a))
end

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(a::ResElem)
   parent(a)(-data(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) + data(b))
end

function -(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) - data(b))
end

function *(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) * data(b))
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

*(a::ResElem, b::Union{Integer, Rational, AbstractFloat}) = parent(a)(data(a) * b)

*(a::ResElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) * b)

*(a::Union{Integer, Rational, AbstractFloat}, b::ResElem) = parent(b)(a * data(b))

*(a::T, b::ResElem{T}) where {T <: RingElem} = parent(b)(a * data(b))

+(a::ResElem, b::Union{Integer, Rational, AbstractFloat}) = parent(a)(data(a) + b)

+(a::ResElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) + b)

+(a::Union{Integer, Rational, AbstractFloat}, b::ResElem) = parent(b)(a + data(b))

+(a::T, b::ResElem{T}) where {T <: RingElem} = parent(b)(a + data(b))

-(a::ResElem, b::Union{Integer, Rational, AbstractFloat}) = parent(a)(data(a) - b)

-(a::ResElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) - b)

-(a::Union{Integer, Rational, AbstractFloat}, b::ResElem) = parent(b)(a - data(b))

-(a::T, b::ResElem{T}) where {T <: RingElem} = parent(b)(a - data(b))

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::ResElem, b::Int)
   parent(a)(powermod(data(a), b, modulus(a)))
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}

Return `true` if $a == b$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   fl = check_parent(a, b, false)
   !fl && return false
   return data(a) == data(b)
end

@doc Markdown.doc"""
    isequal(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}

Return `true` if $a == b$ exactly, otherwise return `false`. This function is
useful in cases where the data of the residues are inexact, e.g. power series
Only if the power series are precisely the same, to the same precision, are
they declared equal by this function.
"""
function isequal(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   fl = check_parent(a, b, false)
   !fl && return false
   return isequal(data(a), data(b))
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(a::ResElem, b::Union{Integer, Rational, AbstractFloat})

Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::ResElem, b::Union{Integer, Rational, AbstractFloat})
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

@doc Markdown.doc"""
    ==(a::Union{Integer, Rational, AbstractFloat}, b::ResElem)

Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::Union{Integer, Rational, AbstractFloat}, b::ResElem)
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

@doc Markdown.doc"""
    ==(a::ResElem{T}, b::T) where {T <: RingElem}

Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::ResElem{T}, b::T) where {T <: RingElem}
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

@doc Markdown.doc"""
    ==(a::T, b::ResElem{T}) where {T <: RingElem}

Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::T, b::ResElem{T}) where {T <: RingElem}
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
    Base.inv(a::ResElem)

Return the inverse of the element $a$ in the residue ring. If an impossible
inverse is encountered, an exception is raised.
"""
function Base.inv(a::ResElem)
   g, ainv = gcdinv(data(a), modulus(a))
   if g != 1
      error("Impossible inverse in inv")
   end
   return parent(a)(ainv)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   fl, q = divides(a, b)
   if !fl
      error("Impossible inverse in divexact")
   end
   return q
end

function divides(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   if iszero(a)
      return true, a
   end
   A = data(a)
   B = data(b)
   R = parent(a)
   m = modulus(R)
   gb = gcd(B, m)
   ub = divexact(B, gb)
   q, r = divrem(A, gb)
   if !iszero(r)
     return false, b
   end
   ub = divexact(B, gb)
   b1 = invmod(ub, divexact(m, gb))
   rs = R(q)*b1
   return true, rs
end

###############################################################################
#
#   GCD
#
###############################################################################

@doc Markdown.doc"""
    gcd(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}

Return a greatest common divisor of $a$ and $b$ if one exists. This is done
by taking the greatest common divisor of the data associated with the
supplied residues and taking its greatest common divisor with the modulus.
"""
function gcd(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(gcd(gcd(data(a), modulus(a)), data(b)))
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::ResElem{T}) where {T <: RingElement}
   a.data = zero!(a.data)
   return a
end

function mul!(c::ResElem{T}, a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   c.data = mod(data(a)*data(b), modulus(a))
   return c
end

function addeq!(c::ResElem{T}, a::ResElem{T}) where {T <: RingElement}
   c.data = mod(data(c) + data(a), modulus(a))
   return c
end

function add!(c::ResElem{T}, a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   c.data = mod(data(a) + data(b), modulus(a))
   return c
end

###############################################################################
#
#   Random functions
#
###############################################################################

RandomExtensions.maketype(R::ResRing, _) = elem_type(R)

# define rand(make(S, v))
function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:ResElem{T},
                                         <:ResRing{T}}}
              ) where {T}
   S, v = sp[][1:end]
   S(rand(rng, v))
end

function RandomExtensions.make(S::ResRing, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1])
   else
      make(S, make(base_ring(S), vs...))
   end
end

rand(rng::AbstractRNG, S::ResRing, v...) = rand(rng, make(S, v...))

rand(S::ResRing, v...) = rand(Random.GLOBAL_RNG, S, v...)

###############################################################################
#
#   ResidueRing constructor
#
###############################################################################

@doc Markdown.doc"""
    ResidueRing(R::Ring, a::RingElement; cached::Bool=true)

Create the residue ring $R/(a)$ where $a$ is an element of the ring $R$. We
require $a \neq 0$. If `cached == true` (the default) then the resulting
residue ring parent object is cached and returned for any subsequent calls
to the constructor with the same base ring $R$ and element $a$. A modulus
of zero is not supported and throws an exception.
"""
function ResidueRing(R::Ring, a::RingElement; cached::Bool = true)
   # Modulus of zero cannot be supported. E.g. A C library could not be expected to
   # do matrices over Z/0 using a Z/nZ type. The former is multiprecision, the latter not.
   iszero(a) && throw(DomainError(a, "Modulus must be nonzero"))
   T = elem_type(R)

   return Generic.ResRing{T}(R(a), cached)
end

function ResidueRing(R::PolyRing, a::RingElement; cached::Bool = true)
   iszero(a) && throw(DomainError(a, "Modulus must be nonzero"))
   !isunit(leading_coefficient(a)) && throw(DomainError(a, "Non-invertible leading coefficient"))
   T = elem_type(R)

   return Generic.ResRing{T}(R(a), cached)
end
