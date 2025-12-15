
###############################################################################

#   Declaration types
#   LocalizedEuclideanRing / LocalizedEuclideanRingElem
#
###############################################################################

# prime might be product of several primes if localized at several primes, those primes are in array primes
@attributes mutable struct LocalizedEuclideanRing{T} <: AbstractAlgebra.Ring
   base_ring::AbstractAlgebra.Ring
   prime::T
   primes::Vector{T}  # in general, not set.
   comp::Bool  # false: den has to be coprime to prime
               # true:  den can ONLY use prime (and powers)

   function LocalizedEuclideanRing{T}(prime::T, primes::Vector{T}, cached::Bool = true, comp::Bool = false) where {T <: RingElem}
      length(primes) == 0 && error("No element to localize at since array of primes is empty")
      if cached && haskey(LocDict, (parent(prime), prime, comp))
         return LocDict[parent(prime), prime, comp]::LocalizedEuclideanRing{T}
      else
         z = new(parent(prime), prime, primes, comp)
         if cached
            LocDict[parent(prime), prime, comp] = z
         end
         return z
      end
   end
   function LocalizedEuclideanRing{T}(prime::T, cached::Bool = true, comp::Bool = false) where {T <: RingElem}
     is_unit(prime) && error("no-point")
     if cached && haskey(LocDict, (parent(prime), prime, comp))
       return LocDict[parent(prime), prime, comp]::LocalizedEuclideanRing{T}
     else
       r = new()
       r.base_ring = parent(prime)
       r.prime = prime
       r.comp = comp
       if cached
          LocDict[parent(prime), prime, comp] = r
       end
       return r
     end
   end
end


const LocDict = Dict{Tuple{AbstractAlgebra.Ring, RingElement, Bool}, AbstractAlgebra.Ring}()

function isin(a, L::LocalizedEuclideanRing{T}) where {T <: RingElem}
  iszero(a) && return true
  L.comp || (!isone(gcd(denominator(a), prime(L))) && return false)
  L.comp && ppio(denominator(a), prime(L))[1] != denominator(a) && return false
  return true
end

mutable struct LocalizedEuclideanRingElem{T} <: AbstractAlgebra.RingElem
   data::FieldElem
   parent::LocalizedEuclideanRing{T}
   function LocalizedEuclideanRingElem{T}(data::FracElem{T}, par::LocalizedEuclideanRing, checked::Bool = true) where {T <: RingElem}
     checked && (isin(data, par) || error("illegal elt"))
     return new{T}(data,par)
   end
end

#=

Let s be an integer (start with localisations in PIDs = Z)
and S = {x | gcd(x, s) = 1} (for s = p, S = Z \setminus p)

  the localsation L = S^-1 R = {a/b | gcd(b, s) = 1}

This is euclidean under N(a/b) = gcd(a, s^infty)

a_1*s_1/b_1 : a_2*s_2/b_2

  divrem(a_1 * s_1 * b_2, s_2) = q, r  => r < s_2
  a_1 s_1 b_2 = q s_2 + r

a_1*s_1/b_1 =   q/(a_2 b_1) *  a_2 s_2/b_2 + r/(b_1 b_2)


This works....


Now the other one:
  S = { s^i : i}
  L = S^-1 R  = { a/b | gcd(b, s^infty) = b}

a_1/s_1 : a_2/s_2   N(a/s) = |a/gcd(a, s^infty)|

a_1s_2 = q a_2 + r
a_1/s_1 = q/s_1 a_2/s_2 + r/(s1 s2)

===========================================
Poly: deg(N(a/b)), rest the same

===========================================
=#

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

elem_type(::Type{LocalizedEuclideanRing{T}}) where {T} = LocalizedEuclideanRingElem{T}

parent_type(::Type{LocalizedEuclideanRingElem{T}}) where {T} = LocalizedEuclideanRing{T}

base_ring_type(::Type{LocalizedEuclideanRing{T}}) where {T} = parent_type(T)

base_ring(L::LocalizedEuclideanRing) = L.base_ring::base_ring_type(L)

parent(a::LocalizedEuclideanRingElem) = a.parent

characteristic(L::LocalizedEuclideanRing) = characteristic(base_ring(L))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

data(a::LocalizedEuclideanRingElem) = a.data

Base.numerator(a::LocalizedEuclideanRingElem{T}, canonicalise::Bool=true) where {T <: RingElement} = numerator(data(a), canonicalise)

Base.denominator(a::LocalizedEuclideanRingElem{T}, canonicalise::Bool=true) where {T <: RingElement} = denominator(data(a), canonicalise)

prime(L::LocalizedEuclideanRing) = L.prime

zero(L::LocalizedEuclideanRing) = L(0)

one(L::LocalizedEuclideanRing) = L(1)

iszero(a::LocalizedEuclideanRingElem) = iszero(data(a))

isone(a::LocalizedEuclideanRingElem) = isone(data(a))

function is_unit(a::LocalizedEuclideanRingElem{T})  where {T <: RingElem}
  return isin(inv(a.data), parent(a))
end

deepcopy_internal(a::LocalizedEuclideanRingElem, dict::IdDict) = parent(a)(deepcopy_internal(data(a), dict))

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, a::LocalizedEuclideanRingElem)
   print(io, data(a))
end

function show(io::IO, L::LocalizedEuclideanRing)
   @show_name(io, L)
   @show_special(io, L)
   io = pretty(io)
   if L.comp
     print(io, "Localization of ", Lowercase(), base_ring(L), " at complement of ", prime(L))
   else
     print(io, "Localization of ", Lowercase(), base_ring(L), " at ", prime(L))
   end
end

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(a::LocalizedEuclideanRingElem)
   parent(a)(-data(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem) where {T}
   check_parent(a,b)
   return LocalizedEuclideanRingElem{T}(data(a) + data(b), parent(a), false)
end

function -(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem{T})  where {T}
   check_parent(a,b)
   return LocalizedEuclideanRingElem{T}(data(a) - data(b), parent(a), false)
end

function *(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem{T})  where {T}
   check_parent(a,b)
   return LocalizedEuclideanRingElem{T}(data(a) * data(b), parent(a), false)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return data(a) == data(b)
end

function Base.hash(a::LocalizedEuclideanRingElem, h::UInt)
   return hash(data(a), h)
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc raw"""
     inv(a::LocalizedEuclideanRingElem{T}, checked::Bool = true)  where {T <: RingElem}
Returns the inverse element of $a$ if $a$ is a unit.
If 'checked = false' the invertibility of $a$ is not checked and the corresponding inverse element
of the Fraction Field is returned.
"""
function Base.inv(a::LocalizedEuclideanRingElem{T}, checked::Bool = true)  where {T}
   b = inv(a.data)
   checked && (isin(b, parent(a)) || error("no unit"))
   return LocalizedEuclideanRingElem{T}(b, parent(a), false)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divides(a::LocalizedEuclideanRingElem, b::LocalizedEuclideanRingElem; checked::Bool = true)
   check_parent(a, b)
   c = divexact(data(a), data(b); check=false)
   isin(c, parent(a)) && return true, parent(a)(c, checked)
   return false, a
end

@doc raw"""
     divexact(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem{T}, checked::Bool = true)  where {T <: RingElem}

Returns element 'c' of given localization s.th. `c`*$b$ = $a$ if such element exists.
If 'checked = true' the result is checked to ensure it is an element of the given
localization.
"""
function divexact(a::LocalizedEuclideanRingElem, b::LocalizedEuclideanRingElem,; checked::Bool=true, check::Bool=true)
   d = divides(a, b; checked=checked)
   d[1] ? d[2] : error("$a not divisible by $b in the given Localization")
end

function Base.divrem(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem{T}, checked::Bool = true)  where {T <: RingElem}
  check_parent(a, b)
  L = parent(a)
  if L.comp
    a1, s1 = ppio(numerator(a.data), L.prime)
    a2, s2 = ppio(numerator(b.data), L.prime)
    b1 = denominator(a)
    b2 = denominator(b)
    q, r = divrem(a1 * s1 * b2, s2)
    return L(q//(a2*b1), checked), L(r//(b1*b2), checked)
  else
    q, r = divrem(numerator(a)*denominator(b), numerator(b))
    return L(q//denominator(a), checked), L(r//(denominator(a)*denominator(b)), checked)
  end
end

function rem(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem{T}) where {T}
  return divrem(a, b)[2]
end

function euclid(a::LocalizedEuclideanRingElem{T}) where {T <: RingElem}
  L = parent(a)
  if L.comp
    return ppio(numerator(a.data), L.prime)[1]
  else
    return ppio(numerator(a.data), L.prime)[2]
  end
end

###############################################################################
#
#   GCD & LCM
#
###############################################################################

function gcd(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem{T}) where {T <: Union{RingElem,Integer}}
   check_parent(a,b)
   iszero(a) && return inv(canonical_unit(b)) * b
   iszero(b) && return inv(canonical_unit(a)) * a
   par = parent(a)
   if par.comp
     elem = ppio(gcd(numerator(a.data), numerator(b.data)), parent(a).prime)[2]
   else
     elem = ppio(gcd(numerator(a.data), numerator(b.data)), parent(a).prime)[1]
   end
   return par(elem)
end

function lcm(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem{T}) where {T <: Union{RingElem,Integer}}
   check_parent(a,b)
   par = parent(a)
   (iszero(a) || iszero(b)) && return par()
   if par.comp
     elem = ppio(lcm(numerator(a.data), numerator(b.data)), parent(a).prime)[2]
   else
     elem = ppio(lcm(numerator(a.data), numerator(b.data)), parent(a).prime)[1]
   end
   return par(elem)
end

###############################################################################
#
#   GCDX
#
###############################################################################

function gcdx(a::LocalizedEuclideanRingElem{T}, b::LocalizedEuclideanRingElem{T}) where {T <: RingElement}
   check_parent(a,b)
   L = parent(a)
   g, u, v = gcdx(numerator(a.data), numerator(b.data))
   c = inv(canonical_unit(L(g)))
   return c*L(g), c*L(u*denominator(a.data)), c*L(v*denominator(b.data))
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::LocalizedEuclideanRingElem, b::Int)
   return parent(a)(data(a)^b)
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{LocalizedEuclideanRingElem{T}}, ::Type{LocalizedEuclideanRingElem{T}}) where {T <: RingElement} = LocalizedEuclideanRingElem{T}

promote_rule(::Type{LocalizedEuclideanRingElem{T}}, ::Type{T}) where {T} = LocalizedEuclideanRingElem{T}

###############################################################################
#
#   Parent object call overload
#
###############################################################################

(L::LocalizedEuclideanRing{T})() where {T} = L(zero(base_ring(L)))

(L::LocalizedEuclideanRing{T})(a::Integer)  where {T} = L(base_ring(L)(a))

function (L::LocalizedEuclideanRing{T})(data::FracElem{T}, checked::Bool = true) where {T <: RingElem}
   return LocalizedEuclideanRingElem{T}(data,L,checked)
end

function (L::LocalizedEuclideanRing{T})(data::Rational{<: Integer}, checked::Bool = true) where {T}
   return LocalizedEuclideanRingElem{T}(base_ring(L)(numerator(data)) // base_ring(L)(denominator(data)),L,checked)
end

function (L::LocalizedEuclideanRing{T})(data::T, checked::Bool = false) where {T <: RingElem}
   return LocalizedEuclideanRingElem{T}(data // parent(data)(1),L,checked)
end

function (L::LocalizedEuclideanRing{T})(a::LocalizedEuclideanRingElem{T}) where {T <: RingElement}
   L != parent(a) && error("No element of $L")
   return a
end

################################################################################
#
#   Valuation
#
################################################################################

@doc raw"""
    valuation(a::LocalizedEuclideanRingElem{T}, p::T) where {T <: RingElement}

Returns the valuation `n` of $a$ at $p$, i.e. the integer `n` s.th $a$ = $p$^`n` * x, where x has valuation 0 at $p$.
"""
valuation(a::LocalizedEuclideanRingElem{T}, p::T) where {T <: RingElement} = valuation(data(a), p)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

@doc raw"""
    canonical_unit(a::LocalizedEuclideanRingElem{T}) where {T <: RingElem}

Returns unit `b`::LocalizedEuclideanRingElem{T} s.th. $a$ * `b` only consists of powers of primes localized at.
"""
function canonical_unit(a::LocalizedEuclideanRingElem{T}) where {T <: RingElem}
   if parent(a).comp
     b = ppio(numerator(a.data), parent(a).prime)[1]
   else
     b = ppio(numerator(a.data), parent(a).prime)[2]
   end
   return parent(a)(b//denominator(a.data))
end

###############################################################################
#
#   Constructors
#
###############################################################################

@doc raw"""
    localization(R::AbstractAlgebra.Ring, prime::T; cached::Bool=true, comp=false) where {T <: RingElement}

Returns the localization of the ring $R$ at the ideal generated by the ring element $prime$. Requires $R$ to
be an euclidean domain and $prime$ to be a prime element, both not checked.
If `cached == true` (the default) then the resulting
localization parent object is cached and returned for any subsequent calls
to the constructor with the same base ring $R$ and element $prime$.
"""
function localization(R::AbstractAlgebra.Ring, prime::T; cached::Bool=true, comp::Bool = false) where {T <: RingElement}
   return LocalizedEuclideanRing{elem_type(R)}(R(prime), cached, comp)
end

@doc raw"""
     localization(R::AbstractAlgebra.Ring, primes::Vector{T}; cached::Bool=true) where {T <: RingElement}

Returns the localization of the ring $R$ at the union of principal ideals that are generated by the ring elements in $primes$.
Requires $R$ to be an euclidean domain and $primes$ to be prime elements, both not checked.
If `cached == true` (the default) then the resulting
localization parent object is cached and returned for any subsequent calls
to the constructor with the same base ring $R$ and elements $primes$.
"""
function localization(R::AbstractAlgebra.Ring, primes::Vector{T}; cached::Bool=true) where {T <: RingElement}
   prime = R(prod(primes))
   return LocalizedEuclideanRing{elem_type(R)}(prime, Vector{elem_type(R)}(primes), cached)
end
