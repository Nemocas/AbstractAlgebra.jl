###############################################################################
#
#   NCRings.jl : Generic not necessarily commutative rings
#
###############################################################################

################################################################################
#
#   Promotion system
#
# The promote_rule functions are not extending Base.promote_rule. The
# AbstractAlgebra promotion system is orthogonal to the built-in julia promotion
# system. The julia system assumes that whenever you have a method signature of
# the form Base.promote_rule(::Type{T}, ::Type{S}) = R, then there is also a
# corresponding Base.convert(::Type{R}, ::T) and similar for S. Since we
# cannot use the julia convert system (we need an instance of the type and not
# the type), we cannot use the julia promotion system.
#
# The AbstractAlgebra promotion system is used to define catch all functions for
# arithmetic between arbitrary ring elements.
#
# TODO: move this to NCRing.jl
#
################################################################################

promote_rule(::Type{T}, ::Type{T}) where T <: NCRingElement = T

function promote_rule_sym(::Type{T}, ::Type{S}) where {T, S}
   U = promote_rule(T, S)
   if U !== Union{}
      return U
   else
      UU = promote_rule(S, T)
      return UU
   end
end

@inline function try_promote(x::S, y::T) where {S <: NCRingElem, T <: NCRingElem}
   U = promote_rule_sym(S, T)
   if S === U
      return true, x, parent(x)(y)
   elseif T === U
      return true, parent(y)(x), y
   else
      return false, x, y
   end
end

function Base.promote(x::S, y::T) where {S <: NCRingElem, T <: NCRingElem}
  fl, u, v = try_promote(x, y)
  if fl
    return u, v
  else
    error("Cannot promote to common type")
  end
end

###############################################################################
#
#   Generic catchall functions
#
###############################################################################

+(x::NCRingElem) = x

+(x::NCRingElem, y::NCRingElem) = +(promote(x, y)...)

+(x::NCRingElem, y::NCRingElement) = x + parent(x)(y)

+(x::NCRingElement, y::NCRingElem) = parent(y)(x) + y

-(x::NCRingElem, y::NCRingElem) = -(promote(x, y)...)

-(x::NCRingElem, y::NCRingElement) = x - parent(x)(y)

-(x::NCRingElement, y::NCRingElem) = parent(y)(x) - y

*(x::NCRingElem) = x

*(x::NCRingElem, y::NCRingElem) = *(promote(x, y)...)

*(x::NCRingElem, y::NCRingElement) = x*parent(x)(y)

*(x::NCRingElement, y::NCRingElem) = parent(y)(x)*y

function ==(x::NCRingElem, y::NCRingElem)
   fl, u, v = try_promote(x, y)
   if fl
     return u == v
   else
     return false
   end
 end
 
==(x::NCRingElem, y::NCRingElement) = x == parent(x)(y)

==(x::NCRingElement, y::NCRingElem) = parent(y)(x) == y

function divexact_left(x::NCRingElem, y::NCRingElem; check::Bool = true)
   return divexact_left(promote(x, y)...)
end

function divexact_right(x::NCRingElem, y::NCRingElem; check::Bool = true)
   return divexact_right(promote(x, y)...)
end

function divexact_left(
   x::NCRingElem,
   y::Union{Integer, Rational, AbstractFloat};
   check::Bool = true)

   return divexact_left(x, parent(x)(y); check = check)
end

function divexact_right(
   x::NCRingElem,
   y::Union{Integer, Rational, AbstractFloat};
   check::Bool = true)

   return divexact_right(x, parent(x)(y); check = check)
end

#Base.:/(x::ModuleElem, y::RingElement) = divexact_right(x, y; check=true)
Base.:/(x::NCRingElem, y::NCRingElement) = divexact_right(x, y; check=true)
Base.:/(x::NCRingElem, y::Union{Integer, Rational, AbstractFloat}) = divexact_right(x, y; check=true)

#Base.:\(y::RingElement, x::ModuleElem) = divexact_left(x, y; check=true)
Base.:\(y::NCRingElement, x::NCRingElem) = divexact_left(x, y; check=true)
Base.:\(y::Union{Integer, Rational, AbstractFloat}, x::NCRingElem) = divexact_left(x, y; check=true)

Base.literal_pow(::typeof(^), x::NCRingElem, ::Val{p}) where {p} = x^p

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc raw"""
    is_unit(a::T) where {T <: NCRingElement}

Return true if $a$ is invertible, else return false.

# Examples
```jldoctest
julia> S, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over rationals, x)

julia> is_unit(x), is_unit(S(1)), is_unit(S(4))
(false, true, true)

julia> is_unit(ZZ(-1)), is_unit(ZZ(4))
(true, false)
```
"""
function is_unit end

@doc raw"""
    is_nilpotent(a::T) where {T <: NCRingElement}

Return true iff $a$ is nilpotent, i.e. a^k == 0 for some k.

# Examples
```jldoctest
julia> R, _ = residue_ring(ZZ,720);

julia> S, x = polynomial_ring(R, :x);

julia> is_nilpotent(30*x), is_nilpotent(30+90*x), is_nilpotent(S(15))
(true, true, false)
```
"""
function is_nilpotent(a::T) where {T <: NCRingElement}
  is_domain_type(T) && return is_zero(a)
  throw(NotImplementedError(:is_nilpotent, a))
end

ConformanceTests._implements(::Type{T}, f::typeof(is_nilpotent)) where {T <: NCRingElement} = is_domain_type(T) || _implements_directly(T, f)


###############################################################################
#
#   Characteristic
#
###############################################################################

function characteristic(R::NCRing)
   error("Characteristic not known")
end

###############################################################################
#
#   One and zero
#
###############################################################################

one(x::NCRingElem) = one(parent(x))

zero(x::NCRingElem) = zero(parent(x))

################################################################################
#
#   Transpose for ring elements
#
################################################################################

transpose(x::T) where {T <: NCRingElem} = deepcopy(x)

###############################################################################
#
#   Delayed reduction
#
###############################################################################

# Fall back to ordinary multiplication
function mul_red!(a::T, b::T, c::T, flag::Bool) where T <: NCRingElement
   return mul!(a, b, c)
end

# Define addmul_delayed_reduction! for all ring elem types
function addmul_delayed_reduction!(a::T, b::T, c::T, d::T) where T <: NCRingElement
   d = mul_red!(d, b, c, false)
   return add!(a, d)
end

# Fall back to nop
function reduce!(a::NCRingElement)
   return a
end

###############################################################################
#
#   Baby-steps giant-steps powering
#
###############################################################################

@doc raw"""
    powers(a::Union{NCRingElement, MatElem}, d::Int)

Return an array $M$ of "powers" of `a` where $M[i + 1] = a^i$ for $i = 0..d$.

# Examples
```jldoctest
julia> M = ZZ[1 2 3; 2 3 4; 4 5 5]
[1   2   3]
[2   3   4]
[4   5   5]

julia> A = powers(M, 4)
5-element Vector{AbstractAlgebra.Generic.MatSpaceElem{BigInt}}:
 [1 0 0; 0 1 0; 0 0 1]
 [1 2 3; 2 3 4; 4 5 5]
 [17 23 26; 24 33 38; 34 48 57]
 [167 233 273; 242 337 394; 358 497 579]
 [1725 2398 2798; 2492 3465 4044; 3668 5102 5957]

```
"""
function powers(a::T, d::Int) where {T <: Union{NCRingElement, MatElem}}
   d < 0 && throw(DomainError(d, "the second argument must be non-negative"))
   a isa MatElem && !is_square(a) && throw(DomainError(a, "matrix must be square"))
   M = Vector{T}(undef, d + 1)
   M[1] = one(a)
   if d > 0
      M[2] = a
      for i = 2:d
         M[i + 1] = M[i] * a
      end
   end
   return M
end

###############################################################################
#
#   oftype
#
################################################################################

# This is a useful fallback when applying methods from Base to our matrices,
# see https://github.com/Nemocas/Nemo.jl/pull/637
Base.oftype(x::NCRingElem, y::NCRingElem) = parent(x)(y)

Base.oftype(x::NCRingElem, y::AbstractFloat) = parent(x)(y)

Base.oftype(x::NCRingElem, y::Integer) = parent(x)(y)

Base.oftype(x::NCRingElem, y::Rational) = parent(x)(y)

###############################################################################
#
#   Scalar for broadcasting
#
###############################################################################

Base.broadcastable(x::NCRingElem) = Ref(x)

################################################################################
#
#   Dot
#
################################################################################

dot(x::NCRingElem, y::NCRingElem) = x * y

dot(x::NCRingElem, y::Union{Integer, Rational, AbstractFloat}) = x * y

dot(x::Union{Integer, Rational, AbstractFloat}, y::NCRingElem) = x * y
