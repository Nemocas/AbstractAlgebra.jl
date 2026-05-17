function ore_extension(R::DT, D::Symbol, δ) where {DT<:Ring}
  σ = sigma_endomorphism(δ)
  OO = ore_extension_only(R,D,δ)
  return OO,gen(OO)
end
ore_extension_only(R::T, D::Symbol, δ) where T<:Ring = OreAlgebra{elem_type(R)}(R, D, δ)

function (R::OreAlgebra{T})() where T
  return zero(R)
end

function (R::OreAlgebra{T})(c::T) where T
  iszero(c) && return zero(R)
  return OreOperator{T}(R,[c],1)
end

function (R::OreAlgebra{T})(c::Union{Integer,Rational,AbstractFloat}) where T
  C = base_ring(R)
  return R(C(c))
end

function (R::OreAlgebra{T})(a::OreOperator{T}) where T
 parent(a) != R && error("Unable to coerce polynomial")
 return a
end

function (R::OreAlgebra{T})(c::Vector{T}) where T<:RingElem
  return OreOperator{T}(R,c,length(c))
end

function (R::OreAlgebra{T})(c::Vector{U}) where {T<:RingElem,U<:RingElem}
  return OreOperator{T}(R,c,length(c))
end

function (R::OreAlgebra{T})(c::Vector{U}) where {T<:RingElem,U<:Integer}
  C = base_ring(R)
  return OreOperator{T}(R,C.(c),length(c))
end

###############################################################################
#
#   OreAlgebra
#
###############################################################################

elem_type(::Type{OreAlgebra{T}}) where T = OreOperator{T}
base_ring_type(::Type{OreAlgebra{T}}) where T = parent_type(T)

base_ring(R::OreAlgebra) = R.base_ring

var(R::OreAlgebra) = R.D
symbols(R::OreAlgebra) = [var(R)]

function gen(R::OreAlgebra{T}) where T
  C = base_ring(R)
  return OreOperator{T}(R,[zero(C),one(C)],2)
end

ngens(R::OreAlgebra) = 1
gens(R::OreAlgebra) = [gen(R)]

function zero(R::OreAlgebra{T}) where T
  return OreOperator{T}(R,T[],0)
end

function one(R::OreAlgebra{T}) where T
  return OreOperator{T}(R,T[one(base_ring(R))],1)
end

derivation(R::OreAlgebra) = R.δ
sigma_endomorphism(R::OreAlgebra) = sigma_endomorphism(derivation(R))

###############################################################################
#
#   OreOperator
#
###############################################################################

parent_type(::Type{OreOperator{T}}) where T = OreAlgebra{T}
parent(a::OreOperator{T}) where T = a.parent

is_domain_type(::Type{OreOperator{T}}) where T = false
is_exact_type(::Type{OreOperator{T}}) where T = is_exact_type(T)

promote_rule(::Type{OreOperator{T}},::Type{OreOperator{T}}) where T<:RingElement = OreOperator{T}
function promote_rule(::Type{OreOperator{T}},::Type{U}) where {T<:RingElement,U<:RingElement}
  promote_rule(T, U) == T ? OreOperator{T} : Union{}
end

function deepcopy_internal(a::OreOperator{T}, dict::IdDict) where T
  C = [deepcopy_internal(c,dict) for c in coefficients(a)]

  return OreOperator{T}(parent(a), C, length(C))
end

coefficients(a::OreOperator{T}) where T = a.coeffs[1:length(a)]
length(a::OreOperator{T}) where T = a.length

iszero(a::OreOperator{T}) where T = a.length == 0
function isone(a::OreOperator{T}) where T
  return length(a) == 1 && first(coefficients(a)) |> isone
end

function canonical_unit(a::OreOperator{T}) where T
  return one(parent(a))
end

function coeff(a::OreOperator{T},i) where T
  if 0 <= i && i <= order(a)
    return a.coeffs[i+1]
  else
    return zero(base_ring(a))
  end
end

setcoeff!(a::OreOperator{T}, n::Int, c::Integer) where T = setcoeff!(a,n,base_ring(a)(c))
function setcoeff!(a::OreOperator{T}, n::Int, c::T) where T
  i = n+1
  if i > length(a.coeffs)
    resize!(a.coeffs,i)
    a.coeffs[a.length+1:i] .= zero(base_ring(a))
  end
  a.coeffs[i] = c
  a.length = normalise(a,length(a.coeffs))
  return a
end

function leading_coefficient(a::OreOperator{T}) where T
  iszero(a) && return base_ring(a)()
  return coefficients(a) |> last
end
function leading_term(a::OreOperator{T}) where T
  iszero(a) && return parent(a)()
  A = parent(a)
  C = base_ring(a)
  c = leading_coefficient(a)
  k = order(a)

  return A([fill(zero(C),k); c])
end
function leading_monomial(a::OreOperator{T}) where T
  iszero(a) && return parent(a)()
  A = parent(a)
  C = base_ring(a)
  k = order(a)

  return A(T[(zero(C) for _ in 1:k)..., one(C)])
end

order(a::OreOperator{T}) where T = iszero(a) ? -1 : length(a)-1

function +(a::OreOperator{T},b::OreOperator{T}) where T<:RingElem
  check_parent(a,b)

  la = length(a)
  lb = length(b)
  l = min(la,lb)
  L = max(la,lb)

  new_c = Vector{T}(undef,L)

  new_c[1:l] .= coefficients(a)[1:l] + coefficients(b)[1:l]
  c = la > lb ? coefficients(a) : coefficients(b)
  new_c[(l+1):L] = c[(l+1):L]

  if all(iszero,new_c)
    return zero(parent(a))
  else
    i = findlast(!iszero,new_c)
    return OreOperator{T}(parent(a), new_c, i)
  end
end

function -(a::OreOperator{T}) where T
  R = parent(a)
  c = -a.coeffs

  return OreOperator{T}(R,c,length(c))
end

function -(a::OreOperator{T},b::OreOperator{T}) where T
  return a + -b
end

function *(a::OreOperator{T},b::T) where T<:RingElem
  return a*parent(a)(b)
end

function *(a::T,b::OreOperator{T}) where T<:RingElem
  R = parent(b)
  res = a.*coefficients(b)
  return OreOperator{T}(R,res,length(b))
end

function *(a::OreOperator{T},b::OreOperator{T}) where T
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
    push!(res,zero(R|>base_ring))
    q = [δ.(q); zero(k)] .+ [zero(k); σ.(q)]
  end

  i = findlast(!iszero,res)
  if isnothing(i)
    return zero(R)
  end

  return OreOperator{T}(R,res[1:i],i)
end

function ^(a::OreOperator{T},i::Int) where T
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


function divexact_right(a::OreOperator{T},b::OreOperator{T}; check=true) where T
  iszero(b) && throw(DivideError())
  q,r = rdivrem(a,b)

  !check && return q
  @req iszero(r) "not an exact division"
  return q
end

function ==(a::OreOperator{T},b::OreOperator{T}) where T
  fl = check_parent(a,b,false)
  !fl && return false
  if a.length != b.length
    return false
  end
  return all(splat(==),zip(a.coeffs,b.coeffs))
end

function ==(a::OreOperator{T}, c::T) where T
  if iszero(c)
    return a.length == 0
  elseif a.length == 1
    return leading_coefficient(a) == c
  else
    return false
  end
end
==(c::T,a::OreOperator{T}) where T = a == c

function is_unit(a::OreOperator{T}) where T
  iszero(a) && return false
  isone(a) && return true
  order(a) == 0 && leading_coefficient(a)|>is_unit && return true
end

function zero!(a::OreOperator{T}) where T
  a.coeffs = elem_type(base_ring(a))[]
  a.length = 0

  return a
end

function one!(a::OreOperator{T}) where T
  a.coeffs = [one(base_ring(a))]
  a.length = 1

  return a
end

function set_length!(a::OreOperator{T},n::Int) where T
  resize!(a.coeffs, n)
  a.length = n

  return a
end

function normalise(a::OreOperator{T}, n::Int) where T
  n = min(n,length(a.coeffs))
  i = findlast(!iszero, a.coeffs[1:n])
  isnothing(i) && return 0
  return i
end

function fit!(a::OreOperator{T}, n::Int) where T
  if n > length(a)
    old_n = length(a)+1
    resize!(a.coeffs, n)
    a.coeffs[old_n:n] .= base_ring(a)|>zero
  end

  return
end


