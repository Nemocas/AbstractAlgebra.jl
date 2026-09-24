###############################################################################
#
#   OrePolyRing
#
###############################################################################
#
#   Constructors
#
###############################################################################

function ore_extension(R::DT, D::Symbol, δ::S) where {DT<:Ring,S<:Map(SkewDerivation)}
  σ = sigma_endomorphism(δ)
  OO = ore_extension_only(R,D,δ)
  return OO,gen(OO)
end
ore_extension_only(R::T, D::Symbol, δ::S) where {T<:Ring,S<:Map(SkewDerivation)} = OrePolyRing{elem_type(R),S}(R, D, δ)

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

elem_type(::Type{OrePolyRing{T,S}}) where {T,S} = OrePolyRingElem{T,S}
base_ring_type(::Type{OrePolyRing{T,S}}) where {T,S} = parent_type(T)
base_ring(R::OrePolyRing) = R.base_ring

var(R::OrePolyRing) = R.D
symbols(R::OrePolyRing) = [var(R)]

function gen(R::OrePolyRing{T,S}) where {T,S}
  C = base_ring(R)
  return OrePolyRingElem{T,S}(R,[zero(C),one(C)],2)
end

AbstractAlgebra.ngens(R::OrePolyRing) = 1
gens(R::OrePolyRing) = [gen(R)]

###############################################################################
#
#   Basic manipulation of rings
#
###############################################################################

function zero(R::OrePolyRing{T,S}) where {T,S}
  return OrePolyRingElem{T,S}(R,T[],0)
end

function one(R::OrePolyRing{T,S}) where {T,S}
  return OrePolyRingElem{T,S}(R,T[one(base_ring(R))],1)
end

derivation(R::OrePolyRing) = R.d
sigma_endomorphism(R::OrePolyRing) = sigma_endomorphism(derivation(R))

###############################################################################
#
#   OrePolyRingElem
#
###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{OrePolyRingElem{T,S}}) where {T,S} = OrePolyRing{T,S}
parent(a::OrePolyRingElem{T,S}) where {T,S} = a.parent

promote_rule(::Type{OrePolyRingElem{T,S}},::Type{OrePolyRingElem{T,S}}) where {T<:RingElement,S} = OrePolyRingElem{T,S}
function promote_rule(::Type{OrePolyRingElem{T,S}},::Type{U}) where {T<:RingElement,U<:RingElement,S}
  promote_rule(T, U) == T ? OrePolyRingElem{T,S} : Union{}
end

function deepcopy_internal(a::OrePolyRingElem{T,S}, dict::IdDict) where {T,S}
  C = [deepcopy_internal(c,dict) for c in coefficients(a)]

  return OrePolyRingElem{T,S}(parent(a), C, length(C))
end

###############################################################################
#
#   Constructors
#
###############################################################################

function (R::OrePolyRing)()
  return zero(R)
end

function (R::OrePolyRing{T,S})(c::T) where {T,S}
  iszero(c) && return zero(R)
  return OrePolyRingElem{T,S}(R,[c],1)
end

function (R::OrePolyRing)(c::Union{Integer,Rational,AbstractFloat})
  C = base_ring(R)
  return R(C(c))
end

function (R::OrePolyRing{T,S})(a::OrePolyRingElem{T,S}) where {T,S}
 parent(a) != R && error("Unable to coerce polynomial")
 return a
end

function (R::OrePolyRing{T,S})(c::Vector{T}) where {T<:RingElem,S}
  return OrePolyRingElem{T,S}(R,c,length(c))
end

function (R::OrePolyRing{T,S})(c::Vector{U}) where {T<:RingElem,U<:RingElem,S}
  return OrePolyRingElem{T,S}(R,c,length(c))
end

function (R::OrePolyRing{T,S})(c::Vector{U}) where {T<:RingElem,U<:Integer,S}
  C = base_ring(R)
  return OrePolyRingElem{T,S}(R,C.(c),length(c))
end

###############################################################################
#
#   Basic manipulation of elements
#
###############################################################################

iszero(a::OrePolyRingElem) = a.length == 0
function isone(a::OrePolyRingElem)
  return length(a) == 1 && first(coefficients(a)) |> isone
end

canonical_unit(a::OrePolyRingElem) = one(parent(a))

function is_unit(a::OrePolyRingElem)
  iszero(a) && return false
  isone(a) && return true
  order(a) == 0 && leading_coefficient(a)|>is_unit && return true
end

order(a::OrePolyRingElem) = iszero(a) ? -1 : length(a)-1
length(a::OrePolyRingElem) = a.length

function set_length!(a::OrePolyRingElem,n::Int)
  resize!(a.coeffs, n)
  a.length = n

  return a
end

coefficients(a::OrePolyRingElem) = a.coeffs[1:length(a)]

function coeff(a::OrePolyRingElem, i::Int)
  if 0 <= i && i <= order(a)
    return a.coeffs[i+1]
  else
    return zero(base_ring(a))
  end
end

setcoeff!(a::OrePolyRingElem, n::Int, c::Integer) = setcoeff!(a,n,base_ring(a)(c))
function setcoeff!(a::OrePolyRingElem{T,S}, n::Int, c::T) where {T,S}
  i = n+1
  if i > length(a.coeffs)
    resize!(a.coeffs,i)
    a.coeffs[a.length+1:i] .= zero(base_ring(a))
  end
  a.coeffs[i] = c
  a.length = normalise(a,length(a.coeffs))
  return a
end

function normalise(a::OrePolyRingElem, n::Int)
  n = min(n,length(a.coeffs))
  i = findlast(!iszero, a.coeffs[1:n])
  isnothing(i) && return 0
  return i
end

function fit!(a::OrePolyRingElem, n::Int)
  if n > length(a)
    old_n = length(a)+1
    resize!(a.coeffs, n)
    a.coeffs[old_n:n] .= base_ring(a)|>zero
  end

  return
end

function leading_coefficient(a::OrePolyRingElem)
  iszero(a) && return base_ring(a)()
  return coefficients(a) |> last
end
function leading_term(a::OrePolyRingElem)
  iszero(a) && return parent(a)()
  return leading_coefficient(a)*leading_monomial(a)
end
function leading_monomial(a::OrePolyRingElem)
  iszero(a) && return parent(a)()
  A = parent(a)
  C = base_ring(a)
  k = order(a)

  return A([zeros(C,k); one(C)])
end

###############################################################################
#
#   Operations
#
###############################################################################

function -(a::OrePolyRingElem)
  R = parent(a)
  c = -a.coeffs

  return R(c)
end

function +(a::OrePolyRingElem{T,S},b::OrePolyRingElem{T,S}) where {T<:RingElem,S}
  check_parent(a,b)
  
  z = parent(a)()

  L = max(length(a),length(b))
  fit!(z,L)

  if L == length(b)
    add!(z,b,a)
  else
    add!(z,a,b)
  end

  n = normalise(z,L)

  if iszero(n) 
    return zero(parent(a))
  else
    return set_length!(z,n)
  end
end

function add!(z::OrePolyRingElem{T},a::OrePolyRingElem{T},b::OrePolyRingElem{T}) where T
  # For this function we assume that la > lb, and that this fits into z
  la = length(a)
  lb = length(b)

  z.coeffs[1:lb] .= coefficients(a)[1:lb] + coefficients(b)
  z.coeffs[(lb+1):la] .= coefficients(a)[(lb+1):la]
  z.length = la

  return z
end

function +(a::OrePolyRingElem,c::Integer)
  C = coefficients(a)
  C[1] += c

  return parent(a)(C)
end
+(c::Integer,a::OrePolyRingElem) = a + c

-(a::OrePolyRingElem,c::Integer) = a + -c
-(c::Integer,a::OrePolyRingElem) = a -  c

function -(a::OrePolyRingElem{T,S},b::OrePolyRingElem{T,S}) where {T,S}
  return a + -b
end

function *(a::OrePolyRingElem,c::Integer)
  C = coefficients(a)
  return parent(a)(c*C)
end
*(c::Integer,a::OrePolyRingElem) = a * c

# This is a non-commutative ring after all, so we do need two different `*` methods

function *(a::OrePolyRingElem{T},b::T) where T<:RingElem
  return a*parent(a)(b)
end

function *(a::T,b::OrePolyRingElem{T,S}) where {T<:RingElem,S}
  R = parent(b)
  res = a.*coefficients(b)
  return OrePolyRingElem{T,S}(R,res,length(b))
end

function _mult_step(k::D,q::Vector{T},d::S) where {D,T,S<:SkewDerivation{D}}
  s = sigma_endomorphism(d)
  return [d.(q); zero(k)] .+ [zero(k); s.(q)]
end
function _mult_step(k::D,q::Vector{T},d::S) where {D,T,S<:TrivialSkewDerivation{D}}
  s = sigma_endomorphism(d)
  return [zero(k); s.(q)]
end

function *(a::OrePolyRingElem{T,S},b::OrePolyRingElem{T,S}) where {T,S}
  check_parent(a,b)
  (iszero(a) || iszero(b)) && return zero(a)

  R = parent(a)
  k = base_ring(R)
  d = derivation(R)

  res = zeros(parent(a)|>base_ring,length(b))
  q = coefficients(b)
  l = length(q)
  for p in coefficients(a)
    res += p.*q
    push!(res,zero(k))
    q = _mult_step(k,q,d)
  end

  i = findlast(!iszero,res)
  if isnothing(i)
    return zero(R)
  end

  return OrePolyRingElem{T,S}(R,res[1:i],i)
end

function ==(a::OrePolyRingElem{T,S},b::OrePolyRingElem{T,S}) where {T,S}
  fl = check_parent(a,b,false)
  !fl && return false
  if a.length != b.length
    return false
  end
  return all(splat(==),zip(a.coeffs,b.coeffs))
end

function ==(a::OrePolyRingElem{T,S}, c::T) where {T,S}
  if iszero(c)
    return a.length == 0
  elseif a.length == 1
    return leading_coefficient(a) == c
  else
    return false
  end
end
==(c::T,a::OrePolyRingElem{T,S}) where {T,S} = a == c

function ^(a::OrePolyRingElem{T,S},i::Int) where {T,S}
  if i < 0
    throw(DomainError(i, "exponent must be >= 0"))
  elseif i == 0
    return one(parent(a))
  elseif i == 1
    return a
  else
    return reduce(*,Iterators.repeated(a,i-1);init=a)
  end
end

function divexact_right(a::OrePolyRingElem{T,S},b::OrePolyRingElem{T,S}; check=true) where {T,S}
  iszero(b) && throw(DivideError())
  q,r = rdivrem(a,b)

  !check && return q
  @req iszero(r) "not an exact division"
  return q
end

###############################################################################
#
#   Unsafe operations
#
###############################################################################

function zero!(a::OrePolyRingElem)
  a.coeffs = elem_type(base_ring(a))[]
  a.length = 0

  return a
end

function one!(a::OrePolyRingElem)
  a.coeffs = [one(base_ring(a))]
  a.length = 1

  return a
end

