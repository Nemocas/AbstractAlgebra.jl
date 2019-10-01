###############################################################################
#
#   Residue.jl : generic residue rings (modulo a principal ideal)
#
###############################################################################

export ResidueRing, inv, modulus, data

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{Res{T}}) where T <: RingElement = ResRing{T}

elem_type(::Type{ResRing{T}}) where {T <: RingElement} = Res{T}

@doc Markdown.doc"""
    base_ring(S::AbstractAlgebra.ResRing{T}) where {T <: RingElement}
> Return the base ring $R$ of the given residue ring $S = R/(a)$.
"""
base_ring(S::AbstractAlgebra.ResRing{T}) where {T <: RingElement} = S.base_ring::parent_type(T)

@doc Markdown.doc"""
    base_ring(r::AbstractAlgebra.ResElem)
> Return the base ring $R$ of the residue ring $R/(a)$ that the supplied
> element $r$ belongs to.
"""
base_ring(r::AbstractAlgebra.ResElem) = base_ring(parent(r))

@doc Markdown.doc"""
    parent(a::AbstractAlgebra.ResElem)
> Return the parent object of the given residue element.
"""
parent(a::AbstractAlgebra.ResElem) = a.parent

isdomain_type(a::Type{T}) where T <: AbstractAlgebra.ResElem = false

function isexact_type(a::Type{T}) where {S <: RingElement, T <: AbstractAlgebra.ResElem{S}}
   return isexact_type(S)
end

function check_parent_type(a::AbstractAlgebra.ResRing{T}, b::AbstractAlgebra.ResRing{T}) where {T <: RingElement}
   # exists only to check types of parents agree
end

function check_parent(a::AbstractAlgebra.ResElem, b::AbstractAlgebra.ResElem, throw::Bool = true)
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

function Base.hash(a::AbstractAlgebra.ResElem, h::UInt)
   b = 0x539c1c8715c1adc2%UInt
   return xor(b, xor(hash(data(a), h), h))
end

@doc Markdown.doc"""
    modulus(R::AbstractAlgebra.ResRing)
> Return the modulus $a$ of the given residue ring $S = R/(a)$.
"""
function modulus(S::AbstractAlgebra.ResRing)
   return S.modulus
end

@doc Markdown.doc"""
    modulus(R::AbstractAlgebra.ResElem)
> Return the modulus $a$ of the residue ring $S = R/(a)$ that the supplied
> residue $r$ belongs to.
"""
function modulus(r::AbstractAlgebra.ResElem)
   return modulus(parent(r))
end

data(a::AbstractAlgebra.ResElem) = a.data

@doc Markdown.doc"""
    zero(R::AbstractAlgebra.ResRing)
> Return the zero element of the given residue ring, i.e. $0 \pmod{a}$ where
> $a$ is the modulus of the residue ring.
"""
zero(R::AbstractAlgebra.ResRing) = R(0)

@doc Markdown.doc"""
    one(R::AbstractAlgebra.ResRing)
> Return $1 \pmod{a}$ where $a$ is the modulus of the residue ring.
"""
one(R::AbstractAlgebra.ResRing) = R(1)

@doc Markdown.doc"""
    iszero(a::AbstractAlgebra.ResElem)
> Return `true` if the supplied element $a$ is zero in the residue ring it
> belongs to, otherwise return `false`.
"""
iszero(a::AbstractAlgebra.ResElem) = iszero(data(a))

@doc Markdown.doc"""
    isone(a::AbstractAlgebra.ResElem)
> Return `true` if the supplied element $a$ is one in the residue ring it
> belongs to, otherwise return `false`.
"""
isone(a::AbstractAlgebra.ResElem) = isone(data(a))

@doc Markdown.doc"""
    isunit(a::AbstractAlgebra.ResElem)
> Return `true` if the supplied element $a$ is invertible in the residue ring
> it belongs to, otherwise return `false`.
"""
function isunit(a::AbstractAlgebra.ResElem)
   g = gcd(data(a), modulus(a))
   return isone(g)
end

deepcopy_internal(a::AbstractAlgebra.ResElem, dict::IdDict) =
   parent(a)(deepcopy(data(a)))

###############################################################################
#
#   Canonicalisation
#
###############################################################################

function canonical_unit(x::AbstractAlgebra.ResElem{<:Union{Integer, RingElem}})
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

function show(io::IO, x::AbstractAlgebra.ResElem)
   print(IOContext(io, :compact => true), data(x))
end

function show(io::IO, a::AbstractAlgebra.ResRing)
   print(IOContext(io, :compact => true), "Residue ring of ", base_ring(a), " modulo ", modulus(a))
end

needs_parentheses(x::AbstractAlgebra.ResElem) = needs_parentheses(data(x))

displayed_with_minus_in_front(x::AbstractAlgebra.ResElem) = displayed_with_minus_in_front(data(x))

show_minus_one(::Type{Res{T}}) where {T <: RingElement} = true

###############################################################################
#
#   Unary operations
#
###############################################################################

@doc Markdown.doc"""
    -(a::AbstractAlgebra.ResElem)
> Return $-a$.
"""
function -(a::AbstractAlgebra.ResElem)
   parent(a)(-data(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

@doc Markdown.doc"""
    +(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
> Return $a + b$.
"""
function +(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) + data(b))
end

@doc Markdown.doc"""
    -(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
> Return $a - b$.
"""
function -(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) - data(b))
end

@doc Markdown.doc"""
    *(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
> Return $a\times b$.
"""
function *(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) * data(b))
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

@doc Markdown.doc"""
    *(a::AbstractAlgebra.ResElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a\times b$.
"""
*(a::AbstractAlgebra.ResElem, b::Union{Integer, Rational, AbstractFloat}) = parent(a)(data(a) * b)

@doc Markdown.doc"""
    *(a::AbstractAlgebra.ResElem{T}, b::T) where {T <: RingElem}
> Return $a\times b$.
"""
*(a::AbstractAlgebra.ResElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) * b)

@doc Markdown.doc"""
    *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.ResElem)
> Return $a\times b$.
"""
*(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.ResElem) = parent(b)(a * data(b))

@doc Markdown.doc"""
    *(a::T, b::AbstractAlgebra.ResElem{T}) where {T <: RingElem}
> Return $a\times b$.
"""
*(a::T, b::AbstractAlgebra.ResElem{T}) where {T <: RingElem} = parent(b)(a * data(b))

@doc Markdown.doc"""
    +(a::AbstractAlgebra.ResElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a + b$.
"""
+(a::AbstractAlgebra.ResElem, b::Union{Integer, Rational, AbstractFloat}) = parent(a)(data(a) + b)

@doc Markdown.doc"""
    +(a::AbstractAlgebra.ResElem{T}, b::T) where {T <: RingElem}
> Return $a + b$.
"""
+(a::AbstractAlgebra.ResElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) + b)

@doc Markdown.doc"""
    +(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.ResElem)
> Return $a + b$.
"""
+(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.ResElem) = parent(b)(a + data(b))

@doc Markdown.doc"""
    +(a::T, b::AbstractAlgebra.ResElem{T}) where {T <: RingElem}
> Return $a + b$.
"""
+(a::T, b::AbstractAlgebra.ResElem{T}) where {T <: RingElem} = parent(b)(a + data(b))

@doc Markdown.doc"""
    -(a::AbstractAlgebra.ResElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a - b$.
"""
-(a::AbstractAlgebra.ResElem, b::Union{Integer, Rational, AbstractFloat}) = parent(a)(data(a) - b)

@doc Markdown.doc"""
    -(a::AbstractAlgebra.ResElem{T}, b::T) where {T <: RingElem}
> Return $a - b$.
"""
-(a::AbstractAlgebra.ResElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) - b)

@doc Markdown.doc"""
    -(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.ResElem)
> Return $a - b$.
"""
-(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.ResElem) = parent(b)(a - data(b))

@doc Markdown.doc"""
    -(a::T, b::AbstractAlgebra.ResElem{T}) where {T <: RingElem}
> Return $a - b$.
"""
-(a::T, b::AbstractAlgebra.ResElem{T}) where {T <: RingElem} = parent(b)(a - data(b))

###############################################################################
#
#   Powering
#
###############################################################################

@doc Markdown.doc"""
    ^(a::AbstractAlgebra.ResElem, b::Int)
> Return $a^b$.
"""
function ^(a::AbstractAlgebra.ResElem, b::Int)
   parent(a)(powmod(data(a), b, modulus(a)))
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
> Return `true` if $a == b$ arithmetically, otherwise return `false`. Recall
> that power series to different precisions may still be arithmetically
> equal to the minimum of the two precisions.
"""
function ==(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   fl = check_parent(a, b, false)
   !fl && return false
   return data(a) == data(b)
end

@doc Markdown.doc"""
    isequal(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
> Return `true` if $a == b$ exactly, otherwise return `false`. This function is
> useful in cases where the data of the residues are inexact, e.g. power series
> Only if the power series are precisely the same, to the same precision, are
> they declared equal by this function.
"""
function isequal(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
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
    ==(a::AbstractAlgebra.ResElem, b::Union{Integer, Rational, AbstractFloat})
> Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::AbstractAlgebra.ResElem, b::Union{Integer, Rational, AbstractFloat})
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

@doc Markdown.doc"""
    ==(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.ResElem)
> Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.ResElem)
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

@doc Markdown.doc"""
    ==(a::AbstractAlgebra.ResElem{T}, b::T) where {T <: RingElem}
> Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::AbstractAlgebra.ResElem{T}, b::T) where {T <: RingElem}
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

@doc Markdown.doc"""
    ==(a::T, b::AbstractAlgebra.ResElem{T}) where {T <: RingElem}
> Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::T, b::AbstractAlgebra.ResElem{T}) where {T <: RingElem}
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
    inv(a::AbstractAlgebra.ResElem)
> Return the inverse of the element $a$ in the residue ring. If an impossible
> inverse is encountered, an exception is raised.
"""
function inv(a::AbstractAlgebra.ResElem)
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

@doc Markdown.doc"""
    divexact(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
> Return $a/b$ where the quotient is expected to be exact.
"""
function divexact(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   fl, q = divides(a, b)
   if !fl
      error("Impossible inverse in divexact")
   end
   return q
end

function divides(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
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
    gcd(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
> Return a greatest common divisor of $a$ and $b$ if one exists. This is done
> by taking the greatest common divisor of the data associated with the
> supplied residues and taking its greatest common divisor with the modulus.
"""
function gcd(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(gcd(gcd(data(a), modulus(a)), data(b)))
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   a.data = zero!(a.data)
   return a
end

function mul!(c::AbstractAlgebra.ResElem{T}, a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   c.data = mod(data(a)*data(b), modulus(a))
   return c
end

function addeq!(c::AbstractAlgebra.ResElem{T}, a::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   c.data = mod(data(c) + data(a), modulus(a))
   return c
end

function add!(c::AbstractAlgebra.ResElem{T}, a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   c.data = mod(data(a) + data(b), modulus(a))
   return c
end

###############################################################################
#
#   Random functions
#
###############################################################################

function rand(rng::AbstractRNG, S::AbstractAlgebra.ResRing{T}, v...) where {T <: RingElement}
   R = base_ring(S)
   return S(rand(rng, R, v...))
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{Res{T}}, ::Type{Res{T}}) where T <: RingElement = Res{T}

function promote_rule(::Type{Res{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? Res{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::ResRing{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::ResRing{T})() where {T <: RingElement}
   z = Res{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function (a::ResRing{T})(b::Integer) where {T <: RingElement}
   z = Res{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function (a::ResRing{T})(b::T) where {T <: RingElem}
   base_ring(a) != parent(b) && error("Operation on incompatible objects")
   z = Res{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function (a::ResRing{T})(b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   a != parent(b) && error("Operation on incompatible objects")
   return b
end

###############################################################################
#
#   ResidueRing constructor
#
###############################################################################

@doc Markdown.doc"""
    ResidueRing(R::AbstractAlgebra.Ring, a::RingElement; cached::Bool=true){T <: RingElement}
> Create the residue ring $R/(a)$ where $a$ is an element of the ring $R$. We
> require $a \neq 0$. If `cached == true` (the default) then the resulting
> residue ring parent object is cached and returned for any subsequent calls
> to the constructor with the same base ring $R$ and element $a$.
"""
function ResidueRing(R::AbstractAlgebra.Ring, a::RingElement; cached::Bool = true)
   iszero(a) && throw(DivideError())
   T = elem_type(R)

   return ResRing{T}(R(a), cached)
end
