###############################################################################
#
#   UniversalRing.jl : Universal ring (variables can be added)
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{UniversalRing{T, U}}) where {T, U} = parent_type(T)
base_ring(S::UniversalRing) = S.base_ring::base_ring_type(S)

coefficient_ring_type(::Type{UniversalRing{T, U}}) where {T, U} = parent_type(U)
coefficient_ring(S::UniversalRing) = coefficient_ring(base_ring(S))::coefficient_ring_type(S)

is_domain_type(::Type{<:UniversalRingElem{T}}) where T = is_domain_type(T)

is_exact_type(::Type{<:UniversalRingElem{T}}) where T = is_exact_type(T)

parent(p::UniversalRingElem) = p.parent::parent_type(p)

elem_type(::Type{UniversalRing{T, U}}) where {T, U} = UniversalRingElem{T, U}

parent_type(::Type{UniversalRingElem{T, U}}) where {T, U} = UniversalRing{T, U}

number_of_variables(S::UniversalRing) = number_of_variables(base_ring(S))

number_of_generators(S::UniversalRing) = number_of_generators(base_ring(S))

symbols(S::UniversalRing) = symbols(base_ring(S))

data(p::UniversalRingElem) = p.p

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(R::UniversalRing{T, U}) where {T, U} = UniversalRingElem{T, U}(zero(base_ring(R)), R)

one(R::UniversalRing{T, U}) where {T, U} = UniversalRingElem{T, U}(one(base_ring(R)), R)

iszero(p::UniversalRingElem) = iszero(data(p))

isone(p::UniversalRingElem) = isone(data(p))

is_unit(p::UniversalRingElem) = is_unit(data(p))

is_zero_divisor(p::UniversalRingElem) = is_zero_divisor(data(p))

function is_zero_divisor_with_annihilator(p::UniversalRingElem{T, U}) where {T, U}
   f, b = is_zero_divisor_with_annihilator(data(p))
   return f, UniversalRingElem{T, U}(b, parent(p))
end

is_gen(p::UniversalRingElem) = is_gen(data(p))

function _ensure_variables(S::UniversalRing, v::Vector{<:VarName})
   idx = Int[]
   current_symbols = symbols(S)
   n = length(current_symbols)
   added_symbols = Symbol[]
   for s_ in v
      s = Symbol(s_)
      i = findfirst(==(s), current_symbols)
      if i === nothing
         push!(added_symbols, s)
         push!(idx, n+length(added_symbols))
      else
         push!(idx, i)
      end
   end
   if !isempty(added_symbols)
      S.base_ring = AbstractAlgebra._add_gens(base_ring(S), added_symbols)
   end
   return idx
end

function gen(S::UniversalRing, s::VarName)
   i = findfirst(==(Symbol(s)), symbols(S))
   if i === nothing
      S.base_ring = AbstractAlgebra._add_gens(base_ring(S), [Symbol(s)])
      i = length(symbols(S))
   end
   return @inbounds gen(S, i)
end

function gen(S::UniversalRing{T, U}, i::Int) where {T, U}
   @boundscheck 1 <= i <= nvars(S) || throw(ArgumentError("generator index out of range"))
   return UniversalRingElem{T, U}(gen(base_ring(S), i), S)
end

function gens(S::UniversalRing{T, U}) where {T, U}
   return [UniversalRingElem{T, U}(g, S) for g in gens(base_ring(S))]
end

# HACK: we abuse the @varnames_interface macro to teach gens for UniversalRing
# some superpowers
function _univ_gens(S::UniversalRing{T, U}, vars::Vector{Symbol}) where {T, U}
   idx = _ensure_variables(S, vars)
   # TRICK: @varnames_interface expects two return values, but we only care
   # for the second; so just return literally nothing for the first
   return nothing, [UniversalRingElem{T, U}(gen(base_ring(S), i), S) for i in idx]
end

AbstractAlgebra.@varnames_interface _univ_gens(R::UniversalRing{T, U}, s) where {T, U}

function gens(S::UniversalRing, a::AbstractAlgebra.VarNames, as::AbstractAlgebra.VarNames...)
   res = _univ_gens(S, a, as...)
   length(res) == 2 && return res[2] # special case for improved backwards compatibility
   return res[2:end]
end

canonical_unit(p::UniversalRingElem) = canonical_unit(data(p))

characteristic(R::UniversalRing) = characteristic(base_ring(R))
is_known(::typeof(characteristic), R::UniversalRing) = is_known(characteristic, base_ring(R))

function Base.hash(p::UniversalRingElem, h::UInt)  # TODO
   b = 0xcf418d4529109236%UInt
   return xor(b, hash(data(p), h))
end

function deepcopy_internal(p::UniversalRingElem, dict::IdDict)
   return UniversalRingElem(deepcopy_internal(data(p), dict), parent(p))
end

Base.copy(f::UniversalRingElem) = UniversalRingElem(copy(data(f)), parent(f))

###############################################################################
#
#   String I/O
#
###############################################################################

universal_ring_name(R::UniversalRing) = "ring"

function show(io::IO, R::UniversalRing)
   @show_name(io, R)
   @show_special(io, R)
   print(io, "Universal $(universal_ring_name(R)) over ")
   show(io, base_ring(base_ring(R)))
end

function expressify(a::UniversalRingElem; context = nothing)
   return expressify(data(a), context = context)
end

@enable_all_show_via_expressify UniversalRingElem

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(p::UniversalRingElem)
   return UniversalRingElem(-data(p), parent(p))
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function univ_promote(x::T, y::T) where {T <: UniversalRingElem}
   check_parent(x, y)
   nx = nvars(parent(data(x)))
   ny = nvars(parent(data(y)))
   if nx == ny
      return x, y
   end
   S = parent(x)
   return S(x), S(y)
end

function +(a::T, b::T) where {T <: UniversalRingElem}
   a, b = univ_promote(a, b)
   return T(data(a) + data(b), parent(a))
end

function -(a::T, b::T) where {T <: UniversalRingElem}
   a, b = univ_promote(a, b)
   return T(data(a) - data(b), parent(a))
end

function *(a::T, b::T) where {T <: UniversalRingElem}
   a, b = univ_promote(a, b)
   return T(data(a)*data(b), parent(a))
end

###############################################################################
#
#   Square root
#
###############################################################################

function Base.sqrt(p::T; check::Bool=true) where {T <: UniversalRingElem}
   S = parent(p)
   s = sqrt(data(p); check=check)
   return T(s, S)
end

is_square(p::UniversalRingElem) = is_square(data(p))

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

for op in (:+, :-, :*)
  @eval begin
    function $op(p::T, n::JuliaRingElement) where {T <: UniversalRingElem}
       S = parent(p)
       return T($op(data(p),n), S)
    end

    function $op(p::UniversalRingElem{T, U}, n::U) where {T<:RingElem, U<:RingElem}
       S = parent(p)
       return UniversalRingElem{T, U}($op(data(p),n), S)
    end

    function $op(n::JuliaRingElement, p::T) where {T <: UniversalRingElem}
       S = parent(p)
       return T($op(n,data(p)), S)
    end

    function $op(n::U, p::UniversalRingElem{T, U}) where {T<:RingElem, U<:RingElem}
       S = parent(p)
       return UniversalRingElem{T, U}($op(n,data(p)), S)
    end
  end
end

function divexact(p::T, n::JuliaRingElement; check::Bool=true) where {T <: UniversalRingElem}
   S = parent(p)
   return T(divexact(data(p), n; check=check), S)
end

function divexact(p::UniversalRingElem{T, U}, n::U; check::Bool=true) where {T<:RingElem, U<:RingElem}
   S = parent(p)
   return UniversalRingElem{T, U}(divexact(data(p), n; check=check), S)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::T, b::T) where {T <: UniversalRingElem}
   check_parent(a, b)
   upgrade!(a)
   upgrade!(b)
   return data(a) == data(b)
end

###############################################################################
#
#   Ad hoc comparison functions
#
###############################################################################

==(p::UniversalRingElem, n::JuliaRingElement) = data(p) == n
==(p::UniversalRingElem{T, U}, n::U) where {T<:RingElem, U<:RingElem} = data(p) == n

==(n::JuliaRingElement, p::UniversalRingElem) = data(p) == n
==(n::U, p::UniversalRingElem{T, U}) where {T<:RingElem, U<:RingElem} = data(p) == n

###############################################################################
#
#   Powering
#
###############################################################################

function ^(p::T, b::Int) where {T <: UniversalRingElem}
   S = parent(p)
   return T(data(p)^b, S)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::T, b::T; check::Bool=true) where {T <: UniversalRingElem}
   a, b = univ_promote(a, b)
   return T(divexact(data(a), data(b); check=check), parent(a))
end

function divides(a::T, b::T) where {T <: UniversalRingElem}
   a, b = univ_promote(a, b)
   flag, q = divides(data(a), data(b))
   return flag, T(q, parent(a))
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function Base.div(a::T, b::T) where {T <: UniversalRingElem}
   a, b = univ_promote(a, b)
   return T(div(data(a), data(b)), parent(a))
end

function Base.divrem(a::T, b::T) where {T <: UniversalRingElem}
   a, b = univ_promote(a, b)
   q, r = divrem(data(a), data(b))
   return T(q, parent(a)), T(r, parent(a))
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(a::T, b::T) where {T <: UniversalRingElem}
   a, b = univ_promote(a, b)
   return T(gcd(data(a), data(b)), parent(a))
end

function lcm(a::T, b::T) where {T <: UniversalRingElem}
   a, b = univ_promote(a, b)
   return T(lcm(data(a), data(b)), parent(a))
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::UniversalRing, _...) = elem_type(S)

function RandomExtensions.make(S::UniversalRing, vs...)
   return Make(S, make(base_ring(S), vs...))
end

function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make2{
                 <:RingElement, <:UniversalRing}})
   S, v = sp[][1:end]
   return UniversalRingElem(rand(rng, v), S)
end

function rand(rng::AbstractRNG, S::UniversalRing, v...)
   return rand(rng, make(S, v...))
end

function rand(S::UniversalRing, v...)
   return rand(Random.default_rng(), S, v...)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::UniversalRingElem)
  a.p = zero!(a.p)
  return a
end

function one!(a::UniversalRingElem)
  a.p = one!(a.p)
  return a
end

function neg!(z::T, a::T) where {T <: UniversalRingElem}
  if parent(data(z)) == parent(data(a))
    z.p = neg!(z.p, a.p)
  else
    z.p = -a.p
  end
  return z
end

function add!(a::T, b::T, c::T) where {T <: UniversalRingElem}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = add!(data(a), data(b), data(c))
  else
    a.p = data(b + c)
  end
  return a
end

function add!(a::T, b::T, c::RingElement) where {T <: UniversalRingElem}
  if parent(data(a)) == parent(data(b))
    a.p = add!(data(a), data(b), c)
  else
    a.p = data(b + c)
  end
  return a
end

add!(a::T, b::RingElement, c::T) where {T <: UniversalRingElem} = add!(a, c, b)

function sub!(a::T, b::T, c::T) where {T <: UniversalRingElem}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = sub!(data(a), data(b), data(c))
  else
    a.p = data(b - c)
  end
  return a
end

function sub!(a::T, b::T, c::RingElement) where {T <: UniversalRingElem}
  if parent(data(a)) == parent(data(b))
    a.p = sub!(data(a), data(b), c)
  else
    a.p = data(b - c)
  end
  return a
end

function sub!(a::T, b::RingElement, c::T) where {T <: UniversalRingElem}
  if parent(data(a)) == parent(data(c))
    a.p = sub!(data(a), b, data(c))
  else
    a.p = data(b - c)
  end
  return a
end

function mul!(a::T, b::T, c::T) where {T <: UniversalRingElem}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = mul!(data(a), data(b), data(c))
  else
    a.p = data(b * c)
  end
  return a
end

function mul!(a::T, b::T, c::RingElement) where {T <: UniversalRingElem}
  if parent(data(a)) == parent(data(b))
    a.p = mul!(data(a), data(b), c)
  else
    a.p = data(b * c)
  end
  return a
end

mul!(a::T, b::RingElement, c::T) where {T <: UniversalRingElem} = mul!(a, c, b)

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{T}, ::Type{T}) where {T <: UniversalRingElem} = T

function promote_rule(::Type{UniversalRingElem{T, U}}, ::Type{V}) where {T, U, V <: RingElement}
   promote_rule(T, V) == T ? UniversalRingElem{T, U} : Union{}
end

###############################################################################
#
#   Parent object overload
#
###############################################################################

function upgrade!(p::UniversalRingElem)
   R = base_ring(parent(p))
   if R != parent(p.p)
      p.p = AbstractAlgebra._upgrade(p.p, R)
   end
   return p
end

function (a::UniversalRing)()
   return UniversalRingElem(base_ring(a)(), a)
end

function (a::UniversalRing)(b::RingElement)
   return UniversalRingElem(base_ring(a)(b), a)
end

function (S::UniversalRing{T, U})(p::UniversalRingElem{T, U}) where {T<:RingElem, U<:RingElement}
   parent(p) !== S && error("Unable to coerce")
   return upgrade!(p)
end

###############################################################################
#
#  Factorization
#
###############################################################################

function _wrap_factorization(f::Fac{T}, S::UniversalRing{T, U}) where {T, U}
   res = Fac{UniversalRingElem{T, U}}()
   res.unit = UniversalRingElem{T, U}(f.unit, S)
   for (fact, expo) in f
      mulpow!(res, UniversalRingElem{T, U}(fact, S), expo)
   end
   return res
end

factor_squarefree(f::UniversalRingElem) = _wrap_factorization(factor_squarefree(data(f)), parent(f))

factor(f::UniversalRingElem) = _wrap_factorization(factor(data(f)), parent(f))
