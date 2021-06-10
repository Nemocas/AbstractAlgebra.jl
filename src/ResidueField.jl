###############################################################################
#
#   ResidueField.jl : residue fields (modulo a principal ideal)
#
###############################################################################

export ResidueField

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring(S::ResField{T}) where {T <: RingElement} = S.base_ring::parent_type(T)

base_ring(r::ResFieldElem) = base_ring(parent(r))

parent(a::ResFieldElem) = a.parent

isdomain_type(a::Type{T}) where T <: ResFieldElem = true

function isexact_type(a::Type{T}) where {S <: RingElement, T <: ResFieldElem{S}}
   return isexact_type(S)
end

function check_parent_type(a::ResField{T}, b::ResField{T}) where {T <: RingElement}
   # exists only to check types of parents agree
end

function check_parent(a::ResFieldElem, b::ResFieldElem, throw::Bool = true)
   if parent(a) != parent(b)
      check_parent_type(parent(a), parent(b))
      fl = modulus(parent(a)) != modulus(parent(b))
      fl && throw && error("Incompatible moduli in residue operation") #CF: maybe extend to divisibility?
      return !fl
   end
   return true
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::ResFieldElem, h::UInt)
   b = 0x539c1c8715c1adc2%UInt
   return xor(b, xor(hash(data(a), h), h))
end

@doc Markdown.doc"""
    modulus(S::ResField)

Return the modulus $a$ of the given residue ring $S = R/(a)$.
"""
function modulus(S::ResField)
   return S.modulus
end

@doc Markdown.doc"""
    modulus(r::ResFieldElem)

Return the modulus $a$ of the residue ring $S = R/(a)$ that the supplied
residue $r$ belongs to.
"""
function modulus(r::ResFieldElem)
   return modulus(parent(r))
end

@doc Markdown.doc"""
    characteristic(r::ResField)

Return the modulus $a$ of the residue ring $S = R/(a)$ that the supplied
residue $r$ belongs to.
"""
function characteristic(r::ResField)
   R = base_ring(r)
   return characteristic(R)
end

@doc Markdown.doc"""
    characteristic(r::ResField{T}) where T <: Integer

Return the modulus $a$ of the residue ring $S = R/(a)$ that the supplied
residue $r$ belongs to.
"""
function characteristic(r::ResField{T}) where T <: Integer
   return modulus(r)
end

data(a::ResFieldElem) = a.data

zero(R::ResField) = R(0)

one(R::ResField) = R(1)

iszero(a::ResFieldElem) = iszero(data(a))

isone(a::ResFieldElem) = isone(data(a))

function isunit(a::ResFieldElem)
   g = gcd(data(a), modulus(a))
   return isone(g)
end

deepcopy_internal(a::ResFieldElem, dict::IdDict) =
   parent(a)(deepcopy(data(a)))

###############################################################################
#
#   Canonicalisation
#
###############################################################################

function canonical_unit(x::ResFieldElem{<:Union{Integer, RingElem}})
  if iszero(x)
    return one(parent(x))
  end
  return x
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(@nospecialize(a::ResFieldElem); context = nothing)
   return expressify(data(a), context = context)
end

function show(io::IO, x::ResFieldElem)
   print(IOContext(io, :compact => true), data(x))
end

function show(io::IO, a::ResField)
   print(IOContext(io, :compact => true), "Residue field of ", base_ring(a), " modulo ", modulus(a))
end

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(a::ResFieldElem)
   parent(a)(-data(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) + data(b))
end

function -(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) - data(b))
end

function *(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) * data(b))
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

*(a::ResFieldElem, b::Union{Integer, Rational, AbstractFloat}) = parent(a)(data(a) * b)

*(a::ResFieldElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) * b)

*(a::Union{Integer, Rational, AbstractFloat}, b::ResFieldElem) = parent(b)(a * data(b))

*(a::T, b::ResFieldElem{T}) where {T <: RingElem} = parent(b)(a * data(b))

+(a::ResFieldElem, b::Union{Integer, Rational, AbstractFloat}) = parent(a)(data(a) + b)

+(a::ResFieldElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) + b)

+(a::Union{Integer, Rational, AbstractFloat}, b::ResFieldElem) = parent(b)(a + data(b))

+(a::T, b::ResFieldElem{T}) where {T <: RingElem} = parent(b)(a + data(b))

-(a::ResFieldElem, b::Union{Integer, Rational, AbstractFloat}) = parent(a)(data(a) - b)

-(a::ResFieldElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) - b)

-(a::Union{Integer, Rational, AbstractFloat}, b::ResFieldElem) = parent(b)(a - data(b))

-(a::T, b::ResFieldElem{T}) where {T <: RingElem} = parent(b)(a - data(b))

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::ResFieldElem, b::Integer)
   parent(a)(powermod(data(a), b, modulus(a)))
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}

Return `true` if $a == b$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   fl = check_parent(a, b, false)
   !fl && return false
   return data(a) == data(b)
end

@doc Markdown.doc"""
    isequal(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}

Return `true` if $a == b$ exactly, otherwise return `false`. This function is
useful in cases where the data of the residues are inexact, e.g. power series
Only if the power series are precisely the same, to the same precision, are
they declared equal by this function.
"""
function isequal(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return isequal(data(a), data(b))
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(a::ResFieldElem, b::Union{Integer, Rational, AbstractFloat})

Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::ResFieldElem, b::Union{Integer, Rational, AbstractFloat})
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

@doc Markdown.doc"""
    ==(a::Union{Integer, Rational, AbstractFloat}, b::ResFieldElem)

Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::Union{Integer, Rational, AbstractFloat}, b::ResFieldElem)
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

@doc Markdown.doc"""
    ==(a::ResFieldElem{T}, b::T) where {T <: RingElem}

Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::ResFieldElem{T}, b::T) where {T <: RingElem}
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

@doc Markdown.doc"""
    ==(a::T, b::ResFieldElem{T}) where {T <: RingElem}

Return `true` if $a == b$ arithmetically, otherwise return `false`.
"""
function ==(a::T, b::ResFieldElem{T}) where {T <: RingElem}
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
    inv(a::ResFieldElem)

Return the inverse of the element $a$ in the residue ring. If an impossible
inverse is encountered, an exception is raised.
"""
function Base.inv(a::ResFieldElem)
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

function divexact(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   fl, q = divides(a, b)
   if !fl
      error("Impossible inverse in divexact")
   end
   return q
end

function divides(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   iszero(b) && error("Division by zero in divides")
   return true, a*inv(b)
end

###############################################################################
#
#   GCD
#
###############################################################################

@doc Markdown.doc"""
    gcd(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}

Return a greatest common divisor of $a$ and $b$ if one exists. This is done
by taking the greatest common divisor of the data associated with the
supplied residues and taking its greatest common divisor with the modulus.
"""
function gcd(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(gcd(gcd(data(a), modulus(a)), data(b)))
end

###############################################################################
#
#   Square root
#
###############################################################################

@doc Markdown.doc"""
    issquare(a::ResFieldElem{T}) where T <: Integer

Return `true` if $a$ is a square.
"""
function issquare(a::ResFieldElem{T}) where T <: Integer
   if iszero(a)
      return true
   end
   p = modulus(a)
   pm1div2 = div(p - 1, 2)
   return isone(a^pm1div2)
end

@doc Markdown.doc"""
    sqrt(a::ResFieldElem{T}) where T <: Integer

Return the square root of $a$ if it is a square, otherwise an exception is
raised.
"""
function Base.sqrt(a::ResFieldElem{T}) where T <: Integer
   U = parent(a)
   p = modulus(a)
   if p == 2 # special case, cannot find a quadratic nonresidue mod 2
      return deepcopy(a)
   end
   # Compute Q, S such that p - 1 = Q*2^S
   Q = p - 1
   S = 0 # power of 2 dividing p - 1
   while iseven(Q)
      Q >>= 1
      S += 1
   end
   # find a quadratic nonresidue z mod p
   z = U(rand(1:p - 1))
   while issquare(z)
      z = U(rand(1:p - 1))
   end
   # set up
   M = S
   c = z^Q
   t = a^Q
   R = a^div(Q + 1, 2)
   # main loop
   while true
      if iszero(t)
         return zero(U)
      end
      if isone(t)
         return R
      end
      u = t
      i = 0
      while i < M
         if isone(u)
            break
         end
         u = u^2
         i += 1
      end
      i == M && error("Not a square in sqrt")
      b = c
      for j = 1:M - i - 1
         b = b^2
      end
      M = i
      c = b^2
      t *= c
      R *= b
   end
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::ResFieldElem{T}) where {T <: RingElement}
   a.data = zero!(a.data)
   return a
end

function mul!(c::ResFieldElem{T}, a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   c.data = mod(data(a)*data(b), modulus(a))
   return c
end

function addeq!(c::ResFieldElem{T}, a::ResFieldElem{T}) where {T <: RingElement}
   c.data = mod(data(c) + data(a), modulus(a))
   return c
end

function add!(c::ResFieldElem{T}, a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   c.data = mod(data(a) + data(b), modulus(a))
   return c
end

###############################################################################
#
#   Random functions
#
###############################################################################

RandomExtensions.maketype(R::ResField, _) = elem_type(R)

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:ResFieldElem{T},
                                         <:ResField{T}}}
              ) where {T}
   S, v = sp[][1:end]
   S(rand(rng, v))
end

function RandomExtensions.make(S::ResField, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1])
   else
      make(S, make(base_ring(S), vs...))
   end
end

rand(rng::AbstractRNG, S::ResField, v...) = rand(rng, make(S, v...))

rand(S::ResField, v...) = rand(Random.GLOBAL_RNG, S, v...)

###############################################################################
#
#   ResidueField constructor
#
###############################################################################

@doc Markdown.doc"""
    ResidueField(R::Ring, a::RingElement; cached::Bool = true)

Create the residue ring $R/(a)$ where $a$ is an element of the ring $R$. We
require $a \neq 0$. If `cached == true` (the default) then the resulting
residue ring parent object is cached and returned for any subsequent calls
to the constructor with the same base ring $R$ and element $a$.
"""
function ResidueField(R::Ring, a::RingElement; cached::Bool = true)
   iszero(a) && throw(DivideError())
   T = elem_type(R)

   return Generic.ResField{T}(R(a), cached)
end

