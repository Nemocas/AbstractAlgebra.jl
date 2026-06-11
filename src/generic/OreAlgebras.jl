###############################################################################
#
#   OrePolyRing
#
###############################################################################
#
#   Constructors
#
###############################################################################

function ore_extension(R::DT, D::Symbol, δ) where {DT<:Ring}
  σ = sigma_endomorphism(δ)
  OO = ore_extension_only(R,D,δ)
  return OO,gen(OO)
end
ore_extension_only(R::T, D::Symbol, δ) where T<:Ring = OrePolyRing{elem_type(R)}(R, D, δ)

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

elem_type(::Type{OrePolyRing{T}}) where T = OrePolyRingElem{T}
base_ring_type(::Type{OrePolyRing{T}}) where T = parent_type(T)
base_ring(R::OrePolyRing) = R.base_ring

var(R::OrePolyRing) = R.D
symbols(R::OrePolyRing) = [var(R)]

function gen(R::OrePolyRing{T}) where T
  C = base_ring(R)
  return OrePolyRingElem{T}(R,[zero(C),one(C)],2)
end

AbstractAlgebra.ngens(R::OrePolyRing) = 1
gens(R::OrePolyRing) = [gen(R)]

###############################################################################
#
#   Basic manipulation of rings
#
###############################################################################

function zero(R::OrePolyRing{T}) where T
  return OrePolyRingElem{T}(R,T[],0)
end

function one(R::OrePolyRing{T}) where T
  return OrePolyRingElem{T}(R,T[one(base_ring(R))],1)
end

derivation(R::OrePolyRing) = R.δ
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

parent_type(::Type{OrePolyRingElem{T}}) where T = OrePolyRing{T}
parent(a::OrePolyRingElem{T}) where T = a.parent

promote_rule(::Type{OrePolyRingElem{T}},::Type{OrePolyRingElem{T}}) where T<:RingElement = OrePolyRingElem{T}
function promote_rule(::Type{OrePolyRingElem{T}},::Type{U}) where {T<:RingElement,U<:RingElement}
  promote_rule(T, U) == T ? OrePolyRingElem{T} : Union{}
end

function deepcopy_internal(a::OrePolyRingElem{T}, dict::IdDict) where T
  C = [deepcopy_internal(c,dict) for c in coefficients(a)]

  return OrePolyRingElem{T}(parent(a), C, length(C))
end

###############################################################################
#
#   Constructors
#
###############################################################################

function (R::OrePolyRing{T})() where T
  return zero(R)
end

function (R::OrePolyRing{T})(c::T) where T
  iszero(c) && return zero(R)
  return OrePolyRingElem{T}(R,[c],1)
end

function (R::OrePolyRing{T})(c::Union{Integer,Rational,AbstractFloat}) where T
  C = base_ring(R)
  return R(C(c))
end

function (R::OrePolyRing{T})(a::OrePolyRingElem{T}) where T
 parent(a) != R && error("Unable to coerce polynomial")
 return a
end

function (R::OrePolyRing{T})(c::Vector{T}) where T<:RingElem
  return OrePolyRingElem{T}(R,c,length(c))
end

function (R::OrePolyRing{T})(c::Vector{U}) where {T<:RingElem,U<:RingElem}
  return OrePolyRingElem{T}(R,c,length(c))
end

function (R::OrePolyRing{T})(c::Vector{U}) where {T<:RingElem,U<:Integer}
  C = base_ring(R)
  return OrePolyRingElem{T}(R,C.(c),length(c))
end

###############################################################################
#
#   Basic manipulation of elements
#
###############################################################################

iszero(a::OrePolyRingElem{T}) where T = a.length == 0
function isone(a::OrePolyRingElem{T}) where T
  return length(a) == 1 && first(coefficients(a)) |> isone
end

function canonical_unit(a::OrePolyRingElem{T}) where T
  return one(parent(a))
end

function is_unit(a::OrePolyRingElem{T}) where T
  iszero(a) && return false
  isone(a) && return true
  order(a) == 0 && leading_coefficient(a)|>is_unit && return true
end

order(a::OrePolyRingElem{T}) where T = iszero(a) ? -1 : length(a)-1
length(a::OrePolyRingElem{T}) where T = a.length

function set_length!(a::OrePolyRingElem{T},n::Int) where T
  resize!(a.coeffs, n)
  a.length = n

  return a
end

coefficients(a::OrePolyRingElem{T}) where T = a.coeffs[1:length(a)]

function coeff(a::OrePolyRingElem{T},i) where T
  if 0 <= i && i <= order(a)
    return a.coeffs[i+1]
  else
    return zero(base_ring(a))
  end
end

setcoeff!(a::OrePolyRingElem{T}, n::Int, c::Integer) where T = setcoeff!(a,n,base_ring(a)(c))
function setcoeff!(a::OrePolyRingElem{T}, n::Int, c::T) where T
  i = n+1
  if i > length(a.coeffs)
    resize!(a.coeffs,i)
    a.coeffs[a.length+1:i] .= zero(base_ring(a))
  end
  a.coeffs[i] = c
  a.length = normalise(a,length(a.coeffs))
  return a
end

function normalise(a::OrePolyRingElem{T}, n::Int) where T
  n = min(n,length(a.coeffs))
  i = findlast(!iszero, a.coeffs[1:n])
  isnothing(i) && return 0
  return i
end

function fit!(a::OrePolyRingElem{T}, n::Int) where T
  if n > length(a)
    old_n = length(a)+1
    resize!(a.coeffs, n)
    a.coeffs[old_n:n] .= base_ring(a)|>zero
  end

  return
end

function leading_coefficient(a::OrePolyRingElem{T}) where T
  iszero(a) && return base_ring(a)()
  return coefficients(a) |> last
end
function leading_term(a::OrePolyRingElem{T}) where T
  iszero(a) && return parent(a)()
  return leading_coefficient(a)*leading_monomial(a)
end
function leading_monomial(a::OrePolyRingElem{T}) where T
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

function -(a::OrePolyRingElem{T}) where T
  R = parent(a)
  c = -a.coeffs

  return OrePolyRingElem{T}(R,c,length(c))
end

function +(a::OrePolyRingElem{T},b::OrePolyRingElem{T}) where T<:RingElem
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

function -(a::OrePolyRingElem{T},b::OrePolyRingElem{T}) where T
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

function *(a::T,b::OrePolyRingElem{T}) where T<:RingElem
  R = parent(b)
  res = a.*coefficients(b)
  return OrePolyRingElem{T}(R,res,length(b))
end

function *(a::OrePolyRingElem{T},b::OrePolyRingElem{T}) where T
  check_parent(a,b)
  (iszero(a) || iszero(b)) && return zero(a)

  R = parent(a)
  k = base_ring(R)
  δ = derivation(R)
  σ = sigma_endomorphism(R)

  res = zeros(parent(a)|>base_ring,length(b))
  q = coefficients(b)
  l = length(q)
  for p in coefficients(a)
    res += p.*q
    push!(res,zero(k))
    q = [δ.(q); zero(k)] .+ [zero(k); σ.(q)]
  end

  i = findlast(!iszero,res)
  if isnothing(i)
    return zero(R)
  end

  return OrePolyRingElem{T}(R,res[1:i],i)
end

function ==(a::OrePolyRingElem{T},b::OrePolyRingElem{T}) where T
  fl = check_parent(a,b,false)
  !fl && return false
  if a.length != b.length
    return false
  end
  return all(splat(==),zip(a.coeffs,b.coeffs))
end

function ==(a::OrePolyRingElem{T}, c::T) where T
  if iszero(c)
    return a.length == 0
  elseif a.length == 1
    return leading_coefficient(a) == c
  else
    return false
  end
end
==(c::T,a::OrePolyRingElem{T}) where T = a == c

function ^(a::OrePolyRingElem{T},i::Int) where T
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

function divexact_right(a::OrePolyRingElem{T},b::OrePolyRingElem{T}; check=true) where T
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

function zero!(a::OrePolyRingElem{T}) where T
  a.coeffs = elem_type(base_ring(a))[]
  a.length = 0

  return a
end

function one!(a::OrePolyRingElem{T}) where T
  a.coeffs = [one(base_ring(a))]
  a.length = 1

  return a
end

