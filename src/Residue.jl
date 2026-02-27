###############################################################################
#
#   Residue.jl : residue rings (modulo a principal ideal)
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{<:ResidueRing{T}}) where T <: RingElement = parent_type(T)

base_ring(S::ResidueRing{T}) where {T <: RingElement} = S.base_ring::parent_type(T)

parent(a::ResElem) = a.parent

is_domain_type(a::Type{T}) where T <: ResElem = false

function is_exact_type(a::Type{T}) where {S <: RingElement, T <: ResElem{S}}
   return is_exact_type(S)
end

function check_parent(a::ResElem, b::ResElem, throw::Bool = true)
   Ra = parent(a)
   Rb = parent(b)
   if Ra != Rb
      fl = typeof(Ra) == typeof(Rb) && modulus(Ra) == modulus(Rb)
      !fl && throw && error("Incompatible moduli in residue operation")
      #CF: maybe extend to divisibility?
      return fl
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

@doc raw"""
    modulus(R::ResidueRing)

Return the modulus $a$ of the given residue ring $S = R/(a)$.
"""
function modulus(S::ResidueRing)
   return S.modulus
end

@doc raw"""
    modulus(R::ResElem)

Return the modulus $a$ of the residue ring $S = R/(a)$ that the supplied
residue $r$ belongs to.
"""
function modulus(r::ResElem)
   return modulus(parent(r))
end

is_trivial(S::ResidueRing) = is_unit(modulus(S))

data(a::ResElem) = a.data

lift(a::ResElem) = data(a)

lift(a::ResElem{Int}) = BigInt(data(a))

zero(R::ResidueRing) = R(0)

one(R::ResidueRing) = R(1)

iszero(a::ResElem) = iszero(data(a))

isone(a::ResElem) = isone(data(a)) || a == one(parent(a))

function is_unit(a::ResElem)
   g = gcd(data(a), modulus(a))
   return isone(g)
end

function is_nilpotent(res::ResElem)
  m = modulus(res)
  r = data(res)
  while true
    g = gcd(r, m)
    (g == m) && return true
    is_one(g) && return false
    m = divexact(m, g)
    g = mod(g, m);   r = g^2  # if computation domain is limited precision integer then g = mod(g,m) guarantees that g^2 will not overflow!
  end
end


# currently residue rings are only allowed over domains
# otherwise this function would be more complicated
is_zero_divisor(a::ResElem) = !is_unit(a)

function is_zero_divisor_with_annihilator(a::ResElem)
   g = gcd(data(a), modulus(a))
   b = divexact(modulus(a), g)  # Modulus must be nonzero, so g is nonzero
   return !is_unit(g), parent(a)(b)
end

annihilator(a::ResElem) = is_zero_divisor_with_annihilator(a)[2]

deepcopy_internal(a::ResElem, dict::IdDict) =
   parent(a)(deepcopy_internal(data(a), dict))

characteristic(R::ResidueRing{T}) where T <: Integer = modulus(R)
is_known(::typeof(characteristic), R::ResidueRing{T}) where T <: Integer = true

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

@enable_all_show_via_expressify ResElem

function show(io::IO, a::ResidueRing)
   @show_name(io, a)
   @show_special(io, a)
   if is_terse(io)
     print(io, "Residue ring")
   else
     io = pretty(io)
     print(io, "Residue ring of ",)
     print(terse(io), Lowercase(), base_ring(a))
     print(io, " modulo ", modulus(a))
   end
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

*(a::ResElem, b::JuliaRingElement) = parent(a)(data(a) * b)

*(a::ResElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) * b)

*(a::JuliaRingElement, b::ResElem) = parent(b)(a * data(b))

*(a::T, b::ResElem{T}) where {T <: RingElem} = parent(b)(a * data(b))

+(a::ResElem, b::JuliaRingElement) = parent(a)(data(a) + b)

+(a::ResElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) + b)

+(a::JuliaRingElement, b::ResElem) = parent(b)(a + data(b))

+(a::T, b::ResElem{T}) where {T <: RingElem} = parent(b)(a + data(b))

-(a::ResElem, b::JuliaRingElement) = parent(a)(data(a) - b)

-(a::ResElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) - b)

-(a::JuliaRingElement, b::ResElem) = parent(b)(a - data(b))

-(a::T, b::ResElem{T}) where {T <: RingElem} = parent(b)(a - data(b))

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::ResElem, b::Int)
   if b < 0
      # powermod throws a DivideError when it should throw an NotInvertibleError
      parent(a)(powermod(data(inv(a)), -b, modulus(a)))
   else
      parent(a)(powermod(data(a), b, modulus(a)))
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc raw"""
    ==(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}

Return `true` if $a == b$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return data(a) == data(b)
end

@doc raw"""
    isequal(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}

Return `true` if $a == b$ exactly, otherwise return `false`. This function is
useful in cases where the data of the residues are inexact, e.g. power series
Only if the power series are precisely the same, to the same precision, are
they declared equal by this function.
"""
function isequal(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   parent(a) == parent(b) || return false
   return isequal(data(a), data(b))
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function ==(a::ResElem, b::JuliaRingElement)
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

function ==(a::JuliaRingElement, b::ResElem)
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

function ==(a::ResElem{T}, b::T) where {T <: RingElem}
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

function ==(a::T, b::ResElem{T}) where {T <: RingElem}
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc raw"""
    Base.inv(a::ResElem)

Return the inverse of the element $a$ in the residue ring. If an impossible
inverse is encountered, an exception is raised.
"""
function Base.inv(a::ResElem)
   g, ainv = gcdinv(data(a), modulus(a))
   isone(g) || throw(NotInvertibleError(a))
   return parent(a)(ainv)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::ResElem{T}, b::ResElem{T}; check::Bool=true) where {T <: RingElement}
   check_parent(a, b)
   fl, q = divides(a, b)
   check && !fl && error("Impossible inverse in divexact")
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

@doc raw"""
    gcd(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}

Return a greatest common divisor of $a$ and $b$ if one exists. This is done
by taking the greatest common divisor of the data associated with the
supplied residues and taking its greatest common divisor with the modulus.
"""
function gcd(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(gcd(gcd(data(a), modulus(a)), data(b)))
end

# Return g, s, t, u, v in R with sv - tu a unit and
#   [s t] [a] = [g]
#   [u v] [b]   [0]
# Generic implementation which uses HNF over the base ring.
# g might not coincide with gcd(a, b) because gcd(a, b) is
# gcd(gcd(data(a), modulus(a)), data(b)) and g is just
# gcd(data(a), data(b)).
function gcdxx(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
  check_parent(a, b)
  R = parent(a)
  M = matrix(base_ring(R), 2, 1, [data(a), data(b)])
  H, U = hermite_form_with_transformation(M)
  @assert is_zero(H[2, 1])
  return R(H[1, 1]), R(U[1, 1]), R(U[1, 2]), R(U[2, 1]), R(U[2, 2])
end

# The operation "Quo" on p. 13 of Storjohann "Algorithms for matrix canonical forms"
function _div_for_howell_form(a::ResElem{T}, b::ResElem{T}) where {T <: RingElement}
  check_parent(a, b)
  return parent(a)(div(data(a), data(b)))
end

# Fallback for euclidean rings (that is, rings implementing the euclidean ring
# interface)
_div_for_howell_form(a::T, b::T) where {T <: RingElement} = div(a, b)

###############################################################################
#
#   Random functions
#
###############################################################################

RandomExtensions.maketype(R::ResidueRing, _) = elem_type(R)

# define rand(make(S, v))
function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:ResElem{T},
                                         <:ResidueRing{T}}}
              ) where {T}
   S, v = sp[][1:end]
   S(rand(rng, v))
end

function RandomExtensions.make(S::ResidueRing, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1])
   else
      Make(S, make(base_ring(S), vs...))
   end
end

rand(rng::AbstractRNG, S::ResidueRing, v...) = rand(rng, make(S, v...))

rand(S::ResidueRing, v...) = rand(Random.default_rng(), S, v...)

###############################################################################
#
#   residue_ring constructor
#
###############################################################################

@doc raw"""
    residue_ring(R::Ring, a::RingElement; cached::Bool=true)

Create the residue ring $R/(a)$ where $a$ is an element of the ring $R$. We
require $a \neq 0$. If `cached == true` (the default) then the resulting
residue ring parent object is cached and returned for any subsequent calls
to the constructor with the same base ring $R$ and element $a$. A modulus
of zero is not supported and throws an exception.
"""
function residue_ring(R::Ring, a::RingElement; cached::Bool = true)
   @req !is_trivial(R) "Zero rings are currently not supported as base ring."
   iszero(a) && throw(DomainError(a, "Modulus must be nonzero"))
   T = elem_type(R)
   S = EuclideanRingResidueRing{T}(R(a), cached)
   return S, Generic.EuclideanRingResidueMap(R, S)
end

function residue_ring(R::PolyRing, a::RingElement; cached::Bool = true)
   @req !is_trivial(R) "Zero rings are currently not supported as base ring."
   iszero(a) && throw(DomainError(a, "Modulus must be nonzero"))
   !is_unit(leading_coefficient(a)) && throw(DomainError(a, "Non-invertible leading coefficient"))
   T = elem_type(R)
   S = EuclideanRingResidueRing{T}(R(a), cached)
   return S, Generic.EuclideanRingResidueMap(R, S)
end

@doc raw"""
    quo(R::Ring, a::RingElement; cached::Bool = true)

Returns `S, f` where `S = residue_ring(R, a)` and `f` is the
projection map from `R` to `S`. This map is supplied as a map with section
where the section is the lift of an element of the residue field back
to the ring `R`.
"""
function quo(R::Ring, a::RingElement; cached::Bool = true)
   S, f = residue_ring(R, a; cached = cached)
   return S, f
end

