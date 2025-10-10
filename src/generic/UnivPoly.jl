###############################################################################
#
#   UnivPoly.jl : Universal polynomial ring (variables can be added)
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{<:UniversalPolyRing{T}}) where T = parent_type(T)
base_ring(S::UniversalPolyRing) = base_ring(mpoly_ring(S))::base_ring_type(S)

coefficient_ring_type(T::Type{<:UniversalPolyRing}) = base_ring_type(T)
coefficient_ring(S::UniversalPolyRing) = base_ring(S)

function is_domain_type(::Type{<:UnivPoly{S}}) where {S <: RingElement}
   return is_domain_type(S)
end

function is_exact_type(::Type{<:UnivPoly{S}}) where {S <: RingElement}
   return is_exact_type(S)
end

parent(p::UnivPoly) = p.parent::parent_type(p)

elem_type(::Type{UniversalPolyRing{T}}) where {T<:RingElement} = UnivPoly{T}

parent_type(::Type{UnivPoly{T}}) where {T<:RingElement} = UniversalPolyRing{T}

function mpoly_ring(S::UniversalPolyRing{T}) where {T<:RingElement}
  return S.mpoly_ring::mpoly_ring_type(T)
end

number_of_variables(S::UniversalPolyRing) = number_of_variables(mpoly_ring(S))

number_of_generators(S::UniversalPolyRing) = number_of_generators(mpoly_ring(S))

symbols(S::UniversalPolyRing) = symbols(mpoly_ring(S))

function vars(p::UnivPoly{T}) where {T}
   S = parent(p)
   V = vars(data(p))
   return [UnivPoly{T}(v, S) for v in V]
end

internal_ordering(p::UniversalPolyRing) = internal_ordering(mpoly_ring(p))

data(p::UnivPoly{T}) where {T<:RingElement} = p.p::mpoly_type(T)

###############################################################################
#
#   Manipulating terms and monomials
#
###############################################################################

exponent_vector(p::UnivPoly, i::Int) = exponent_vector(data(p), i)

function exponent(p::UnivPoly, i::Int, j::Int)
   if j > nvars(parent(data(p))) && j <= nvars(parent(p))
      return 0
   end
   return exponent(data(p), i, j)
end

function set_exponent_vector!(p::UnivPoly, i::Int, exps::Vector{Int})
   S = parent(p)
   len = length(exps)
   if len != nvars(parent(data(p)))
      p.p = upgrade(S, data(p))
      if len < nvars(S)
         exps = vcat(exps, zeros(Int, nvars(S) - len))
      end
   end
   p.p = set_exponent_vector!(data(p), i, exps)
   return p
end

function coeff(p::UnivPoly, exps::Vector{Int})
   S = parent(p)
   len = length(exps)
   n = nvars(parent(data(p)))
   if len > n
      if !iszero(exps[n + 1:len])
         return base_ring(S)()
      end
      return coeff(data(p), exps[1:n])
   end
   n = nvars(parent(data(p)))
   if len < n
      exps = vcat(exps, zeros(Int, n - len))
   end
   return coeff(data(p), exps)
end

function setcoeff!(p::UnivPoly, exps::Vector{Int}, c::T) where T <: RingElement
   c = base_ring(data(p))(c)
   S = parent(p)
   len = length(exps)
   if len != nvars(parent(data(p)))
      p.p = upgrade(S, data(p))
      if len < nvars(S)
         exps = vcat(exps, zeros(Int, nvars(S) - len))
      end
   end
   p.p = setcoeff!(data(p), exps, c)
   return p
end

function setcoeff!(a::UnivPoly{T}, i::Int, c::RingElement) where {T}
   setcoeff!(data(a), i, c)
   return a
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(R::UniversalPolyRing{T}) where {T} = UnivPoly{T}(zero(mpoly_ring(R)), R)

one(R::UniversalPolyRing{T}) where {T} = UnivPoly{T}(one(mpoly_ring(R)), R)

iszero(p::UnivPoly) = iszero(data(p))

isone(p::UnivPoly) = isone(data(p))

is_unit(p::UnivPoly) = is_unit(data(p))

is_zero_divisor(p::UnivPoly) = is_zero_divisor(data(p))

function is_zero_divisor_with_annihilator(p::UnivPoly{T}) where {T}
   f, b = is_zero_divisor_with_annihilator(data(p))
   return f, UnivPoly{T}(b, parent(p))
end

is_gen(p::UnivPoly) = is_gen(data(p))

is_homogeneous(p::UnivPoly) = is_homogeneous(data(p))

is_monomial(p::UnivPoly) = is_monomial(data(p))

is_constant(p::UnivPoly) = is_constant(data(p))

is_term(p::UnivPoly) = is_term(data(p))

coeff(p::UnivPoly, i::Int) = coeff(data(p), i)

function coeff(p::UnivPoly{T}, m::UnivPoly{T}) where {T}
   !is_monomial(m) && error("Not a monomial")
   check_parent(p, m)
   v1 = first(exponent_vectors(m))
   len = length(v1)
   n = nvars(parent(data(p)))
   R = base_ring(p)
   if len > n
      if !iszero(v1[n + 1:len])
         return zero(R)
      end
      v1 = v1[1:n]
   end
   if len < n
      v1 = vcat(v1, zeros(Int, n - len))
   end
   cvzip = zip(coefficients(data(p)), exponent_vectors(data(p)))
   for (c, v) in cvzip
      if v == v1
         return c
      end
   end
   return zero(R)
end

function monomial(p::UnivPoly{T}, i::Int) where {T}
   S = parent(p)
   m = monomial(data(p), i)
   return UnivPoly{T}(m, S)
end

function monomial!(m::UnivPoly{T}, p::UnivPoly{T}, i::Int) where {T}
   parent(m) != parent(p) && error("Incompatible monomial")
   if parent(data(m)) != parent(data(p))
      m.p = parent(data(p))()
   end
   m.p = monomial!(data(m), data(p), i)
   return m
end

function term(p::UnivPoly{T}, i::Int) where {T}
   S = parent(p)
   t = term(data(p), i)
   return UnivPoly{T}(t, S)
end

leading_coefficient(p::UnivPoly) = leading_coefficient(data(p))

trailing_coefficient(p::UnivPoly) = trailing_coefficient(data(p))

function tail(p::UnivPoly{T}) where {T}
   S = parent(p)
   return UnivPoly{T}(tail(data(p)), S)
end

constant_coefficient(p::UnivPoly) = constant_coefficient(data(p))

function leading_monomial(p::UnivPoly{T}) where {T}
   S = parent(p)
   return UnivPoly{T}(leading_monomial(data(p)), S)
end

function leading_term(p::UnivPoly{T}) where {T}
   S = parent(p)
   return UnivPoly{T}(leading_term(data(p)), S)
end

max_fields(p::UnivPoly) = max_fields(data(p))

function degree(p::UnivPoly, i::Int)
   if i <= nvars(parent(p)) && i > nvars(parent(data(p)))
      return 0
   end
   return degree(data(p), i)
end

function degree(f::UnivPoly{T}, x::UnivPoly{T}) where {T}
   check_parent(f, x)
   return degree(f, var_index(x))
end

function degrees(p::UnivPoly)
   v = degrees(data(p))
   len = length(v)
   num = nvars(parent(p))
   if len < num
      v = vcat(v, zeros(Int, num - len))
   end
   return v
end

total_degree(p::UnivPoly) = total_degree(data(p))

length(p::UnivPoly) = length(data(p))

function _ensure_variables(S::UniversalPolyRing, v::Vector{<:VarName})
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
      new_symbols = vcat(current_symbols, added_symbols)
      S.mpoly_ring = AbstractAlgebra.polynomial_ring_only(base_ring(S), new_symbols; internal_ordering=internal_ordering(S), cached=false)
   end
   return idx
end

function gen(S::UniversalPolyRing, s::VarName)
   i = findfirst(==(Symbol(s)), symbols(S))
   if i === nothing
      new_symbols = copy(symbols(S))
      push!(new_symbols, Symbol(s))
      i = length(new_symbols)
      S.mpoly_ring = AbstractAlgebra.polynomial_ring_only(base_ring(S), new_symbols; internal_ordering=internal_ordering(S), cached=false)
   end
   return @inbounds gen(S, i)
end

function gen(S::UniversalPolyRing{T}, i::Int) where {T}
   @boundscheck 1 <= i <= nvars(S) || throw(ArgumentError("generator index out of range"))
   return UnivPoly{T}(gen(mpoly_ring(S), i), S)
end

function gens(S::UniversalPolyRing{T}) where {T}
   n = nvars(S)
   return UnivPoly{T}[gen(S, i) for i in 1:n]
end

# HACK: we abuse the @varnames_interface macro to teach gens for UniversalPolyRing
# some superpowers
function _univ_poly_gens(S::UniversalPolyRing{T}, vars::Vector{Symbol}) where {T}
   idx = _ensure_variables(S, vars)
   # TRICK: @varnames_interface expects two return values, but we only care
   # for the second; so just return literally nothing for the first
   return nothing, [UnivPoly{T}(gen(mpoly_ring(S), i), S) for i in idx]
end

AbstractAlgebra.@varnames_interface _univ_poly_gens(R::UniversalPolyRing{T}, s) where {T}

function gens(S::UniversalPolyRing{T}, a::AbstractAlgebra.VarNames, as::AbstractAlgebra.VarNames...) where {T}
   res = _univ_poly_gens(S, a, as...)
   length(res) == 2 && return res[2] # special case for improved backwards compatibility
   return res[2:end]
end

var_index(x::UnivPoly) = var_index(data(x))

function vars(p::UnivPoly{T}) where {T <: RingElement}
   V = vars(data(p))
   S = parent(p)
   return [UnivPoly{T}(v, S) for v in V]
end

canonical_unit(p::UnivPoly) = canonical_unit(data(p))

characteristic(R::UniversalPolyRing) = characteristic(base_ring(R))

function Base.hash(p::UnivPoly, h::UInt)
   b = 0xcf418d4529109236%UInt
   for (c, v) in zip(coefficients(data(p)), exponent_vectors(data(p)))
      l = length(v)
      while l > 0 && iszero(v[l])
         l -= 1
      end
      b = xor(b, xor(Base.hash(v[1:l], h), h))
      b = xor(b, xor(hash(c, h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

function deepcopy_internal(p::UnivPoly{T}, dict::IdDict) where {T}
   return UnivPoly{T}(deepcopy_internal(data(p), dict), parent(p))
end

Base.copy(f::UnivPoly{T}) where {T} = UnivPoly{T}(copy(data(f)), parent(f))

###############################################################################
#
#   Multivariate coefficients
#
###############################################################################

function coeff(p::UnivPoly{T}, vars::Vector{Int}, exps::Vector{Int}) where {T}
   len = length(vars)
   @req len == length(exps) "Number of variables does not match number of exponents"
   S = parent(p)
   n = nvars(S)
   num = nvars(parent(data(p)))
   vars2 = Vector{Int}(undef, 0)
   exps2 = Vector{Int}(undef, 0)
   for i = 1:len
      @req 1 <= vars[i] <= nvars(S) "Variable index not in range"
      if vars[i] <= num
         push!(vars2, vars[i])
         push!(exps2, exps[i])
      elseif exps[i] > 0
         return zero(S)
      end
   end
   return UnivPoly{T}(coeff(data(p), vars2, exps2), S)
end

function coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where {U, T <: UnivPoly{U}}
   varidx = [var_index(x) for x in vars]
   return coeff(a, varidx, exps)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::UniversalPolyRing)
   @show_name(io, R)
   @show_special(io, R)
   print(io, "Universal Polynomial Ring over ")
   show(io, base_ring(R))
end

function expressify(a::UnivPoly, x = symbols(parent(a)); context = nothing)
   return expressify(data(a), x, context = context)
end

@enable_all_show_via_expressify UnivPoly

###############################################################################
#
#   Unary operations
#
###############################################################################
   
function -(p::UnivPoly{T}) where {T}
   return UnivPoly{T}(-data(p), parent(p))
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function univ_promote(x::UnivPoly{T}, y::UnivPoly{T}) where {T <: RingElement}
   check_parent(x, y)
   nx = nvars(parent(data(x)))
   ny = nvars(parent(data(y)))
   if nx == ny
      return x, y
   end
   S = parent(x)
   return S(x), S(y)
end

function +(a::UnivPoly{T}, b::UnivPoly{T}) where {T}
   a, b = univ_promote(a, b)
   return UnivPoly{T}(data(a) + data(b), parent(a))
end

function -(a::UnivPoly{T}, b::UnivPoly{T}) where {T}
   a, b = univ_promote(a, b)
   return UnivPoly{T}(data(a) - data(b), parent(a))
end

function *(a::UnivPoly{T}, b::UnivPoly{T}) where {T}
   a, b = univ_promote(a, b)
   return UnivPoly{T}(data(a)*data(b), parent(a))
end

###############################################################################
#
#   Iterators
#
###############################################################################

function Base.iterate(x::UnivPolyCoeffs)
   if length(x.poly) >= 1
      return coeff(x.poly, 1), 1
   else
      return nothing
   end
end

function Base.iterate(x::UnivPolyCoeffs, state)
   state += 1
   if length(x.poly) >= state
      return coeff(x.poly, state), state
   else
      return nothing
   end
end

function Base.iterate(x::UnivPolyExponentVectors)
   if length(x.poly) >= 1
      return exponent_vector(x.poly, 1), 1
   else
      return nothing
   end
end

function Base.iterate(x::UnivPolyExponentVectors, state)
   state += 1
   if length(x.poly) >= state
      return exponent_vector(x.poly, state), state
   else
      return nothing
   end
end

function Base.iterate(x::UnivPolyTerms)
   if length(x.poly) >= 1
      return term(x.poly, 1), 1
   else
      return nothing
   end
end

function Base.iterate(x::UnivPolyTerms, state)
   state += 1
   if length(x.poly) >= state
      return term(x.poly, state), state
   else
      return nothing
   end
end

function Base.iterate(x::UnivPolyMonomials)
   if length(x.poly) >= 1
      return monomial(x.poly, 1), 1
   else
      return nothing
   end
end

function Base.iterate(x::UnivPolyMonomials, state)
   state += 1
   if length(x.poly) >= state
      return monomial(x.poly, state), state
   else
      return nothing
   end
end

function Base.length(x::Union{UnivPolyCoeffs, UnivPolyExponentVectors, UnivPolyTerms, UnivPolyMonomials})
   return length(x.poly)
end

function Base.eltype(::Type{UnivPolyCoeffs{T}}) where T <: AbstractAlgebra.UniversalPolyRingElem{S} where S <: RingElement
   return S
end

function Base.eltype(::Type{UnivPolyExponentVectors{T}}) where T <: AbstractAlgebra.UniversalPolyRingElem{S} where S <: RingElement
   return Vector{Int}
end

function Base.eltype(::Type{UnivPolyMonomials{T}}) where T <: AbstractAlgebra.UniversalPolyRingElem{S} where S <: RingElement
   return T
end

function Base.eltype(::Type{UnivPolyTerms{T}}) where T <: AbstractAlgebra.UniversalPolyRingElem{S} where S <: RingElement
   return T
end

###############################################################################
#
#   Square root
#
###############################################################################

function Base.sqrt(p::UnivPoly{T}; check::Bool=true) where {T}
   S = parent(p)
   s = sqrt(data(p); check=check)
   return UnivPoly{T}(s, S)
end

# See generic documentation in NCRings.jl
function is_square(p::UnivPoly)
   return is_square(data(p))
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

for op in (:+, :-, :*)
  @eval begin
    function $op(p::UnivPoly{T}, n::Union{Integer, Rational, AbstractFloat}) where {T}
       S = parent(p)
       return UnivPoly{T}($op(data(p),n), S)
    end

    function $op(p::UnivPoly{T}, n::T) where {T <: RingElem}
       S = parent(p)
       return UnivPoly{T}($op(data(p),n), S)
    end

    function $op(n::Union{Integer, Rational, AbstractFloat}, p::UnivPoly{T}) where {T}
       S = parent(p)
       return UnivPoly{T}($op(n,data(p)), S)
    end

    function $op(n::T, p::UnivPoly{T}) where {T <: RingElem}
       S = parent(p)
       return UnivPoly{T}($op(n,data(p)), S)
    end
  end
end

function divexact(p::UnivPoly{T}, n::Union{Integer, Rational, AbstractFloat}; check::Bool=true) where {T}
   S = parent(p)
   return UnivPoly{T}(divexact(data(p), n; check=check), S)
end

function divexact(p::UnivPoly{T}, n::T; check::Bool=true) where {T <: RingElem}
   S = parent(p)
   return UnivPoly{T}(divexact(data(p), n; check=check), S)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::UnivPoly{T}, b::UnivPoly{T}) where {T}
   check_parent(a, b)

   # quick check if underlying parents agree
   if parent(data(a)) === parent(data(b))
     return data(a) == data(b)
   end

   # check for "generators"
   fl1, i1 = AbstractAlgebra._is_gen_with_index(data(a))
   fl2, i2 = AbstractAlgebra._is_gen_with_index(data(b))
   if fl1 && fl2
     return i1 == i2
   elseif fl1 != fl2
     return false
   end

   if length(a) != length(b)
      return false
   end
   n1 = nvars(parent(data(a)))
   n2 = nvars(parent(data(b)))
   if n1 == n2
      for (v1, v2) in zip(exponent_vectors(data(a)), exponent_vectors(data(b)))
         if v1 != v2
            return false
         end
      end
   elseif n1 > n2
      for (v1, v2) in zip(exponent_vectors(data(a)), exponent_vectors(data(b)))
         if v1[1:n2] != v2 || !iszero(v1[n2 + 1:n1])
            return false
         end
      end
   else # n2 > n1
      for (v1, v2) in zip(exponent_vectors(data(a)), exponent_vectors(data(b)))
         if v1 != v2[1:n1] || !iszero(v2[n1 + 1:n2])
            return false
         end
      end
   end
   for (c1, c2) in zip(coefficients(data(a)), coefficients(data(b)))
      if c1 != c2
         return false
      end
   end
   return true
end

function isless(a::UnivPoly{T}, b::UnivPoly{T}) where {T}
   check_parent(a, b)
   if parent(data(a)) === parent(data(b))
      return isless(data(a), data(b))
   end
   S = parent(a)
   num = nvars(S)
   s = data(a)
   t = data(b)
   if nvars(parent(s)) != num
      s = upgrade(S, s)
   end
   if nvars(parent(t)) != num
      t = upgrade(S, t)
   end
   return isless(s, t)
end

###############################################################################
#
#   Ad hoc comparison functions
#
###############################################################################

==(p::UnivPoly, n::Union{Integer, Rational, AbstractFloat}) = data(p) == n

==(n::Union{Integer, Rational, AbstractFloat}, p::UnivPoly) = data(p) == n

==(p::UnivPoly{T}, n::T) where {T <: RingElem} = data(p) == n

==(n::T, p::UnivPoly{T}) where {T <: RingElem} = data(p) == n

###############################################################################
#
#   Powering
#
###############################################################################

function ^(p::UnivPoly{T}, b::Int) where {T}
   S = parent(p)
   return UnivPoly{T}(data(p)^b, S)
end

###############################################################################
#
#   Inflation/deflation
#
###############################################################################

deflation(p::UnivPoly{T}) where {T} = deflation(data(p))

function deflate(p::UnivPoly{T}, shift::Vector{Int}, defl::Vector{Int}) where {T}
   S = parent(p)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   num = nvars(parent(data(p)))
   pp = data(p)
   if vlen == num
      return UnivPoly{T}(deflate(pp, shift, defl), S)
   end
   if vlen > num
      pp = upgrade(S, pp)
      num = nvars(parent(pp))
   end
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return UnivPoly{T}(deflate(pp, shift, defl), S)
end

function deflate(p::UnivPoly{T}, defl::Vector{Int}) where {T}
   return deflate(p, zeros(Int, length(defl)), defl)
end

function inflate(p::UnivPoly{T}, shift::Vector{Int}, defl::Vector{Int}) where {T}
   S = parent(p)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   num = nvars(parent(data(p)))
   pp = data(p)
   if vlen == num
      return UnivPoly{T}(inflate(pp, shift, defl), S)
   end
   if vlen > num
      pp = upgrade(S, pp)
      num = nvars(parent(pp))
   end
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return UnivPoly{T}(inflate(pp, shift, defl), S)
end

function inflate(p::UnivPoly{T}, defl::Vector{Int}) where {T}
   return inflate(p, zeros(Int, length(defl)), defl)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::UnivPoly{T}, b::UnivPoly{T}; check::Bool=true) where {T}
   a, b = univ_promote(a, b)
   return UnivPoly{T}(divexact(data(a), data(b); check=check), parent(a))
end

function divides(a::UnivPoly{T}, b::UnivPoly{T}) where {T}
   a, b = univ_promote(a, b)
   flag, q = divides(data(a), data(b))
   return flag, UnivPoly{T}(q, parent(a))
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function Base.div(a::UnivPoly{T}, b::UnivPoly{T}) where {T}
   a, b = univ_promote(a, b)
   return UnivPoly{T}(div(data(a), data(b)), parent(a))
end

function Base.divrem(a::UnivPoly{T}, b::UnivPoly{T}) where {T}
   a, b = univ_promote(a, b)
   q, r = divrem(data(a), data(b))
   return UnivPoly{T}(q, parent(a)), UnivPoly{T}(r, parent(a))
end

###############################################################################
#
#   Derivative
#
###############################################################################

function derivative(p::UnivPoly{T}, j::Int) where {T}
   j > nvars(parent(p)) && error("No such variable")
   if j > nvars(parent(data(p)))
      return zero(parent(p))
   end
   return UnivPoly{T}(derivative(data(p), j), parent(p))
end

function derivative(p::UnivPoly{T}, x::UnivPoly{T}) where {T}
   return derivative(p, var_index(x))
end

###############################################################################
#
#   Remove and valuation
#
###############################################################################

function remove(z::UnivPoly{T}, p::UnivPoly{T}) where {T}
   z, p = univ_promote(z, p)
   val, q = remove(data(z), data(p))
   return val, UnivPoly{T}(q, parent(z))
end

function valuation(z::UnivPoly{T}, p::UnivPoly{T}) where {T}
  v, _ = remove(z, p)
  return v
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(a::UnivPoly{T}, A::Vector{T}) where {T <: RingElem}
   R = base_ring(a)
   n = length(A)
   num = nvars(parent(data(a)))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
      if nvars(parent(data(a))) == 0
         return constant_coefficient(data(a))*one(parent(A[1]))
      end
      return evaluate(data(a), A[1:num])
   end
   if n < num
      A = vcat(A, [zero(R) for i = 1:num - n])
   end
   return evaluate(data(a), A)
end

function evaluate(a::UnivPoly{T}, A::Vector{V}) where {T <: RingElement, V <: Union{Integer, Rational, AbstractFloat}}
   n = length(A)
   num = nvars(parent(data(a)))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
      if nvars(parent(data(a))) == 0
         return constant_coefficient(data(a))*one(parent(A[1]))
      end
      return evaluate(data(a), A[1:num])
   end
   if n < num
      A = vcat(A, zeros(V, num - n))
   end
   return evaluate(data(a), A)
end

function evaluate(a::UnivPoly{T}, A::Vector{V}) where {T <: RingElement, V <: RingElement}
   n = length(A)
   num = nvars(parent(data(a)))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
      if nvars(parent(data(a))) == 0
         return constant_coefficient(data(a))*one(parent(A[1]))
      end
      return evaluate(data(a), A[1:num])
   end
   if n < num
      if n == 0
         R = base_ring(a)
         return evaluate(data(a), [zero(R) for _ in 1:num])
      else
         R = parent(A[1])
         A = vcat(A, [zero(R) for _ in 1:num-n])
         return evaluate(data(a), A)
     end
   end
   return evaluate(data(a), A)
end

function (a::UnivPoly{T})() where {T <: RingElement}
   return evaluate(a, T[])
end

function (a::UnivPoly{T})(vals::T...) where {T <: RingElement}
   return evaluate(a, [vals...])
end

function (a::UnivPoly{T})(val::V, vals::V...) where {T <: RingElement, V <: Union{Integer, Rational, AbstractFloat}}
   return evaluate(a, [val, vals...])
end

function (a::UnivPoly{T})(vals::Union{NCRingElem, RingElement}...) where {T <: RingElement}
   A = [vals...]
   n = length(vals)
   num = nvars(parent(data(a)))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
      if nvars(parent(data(a))) == 0
         return constant_coefficient(data(a))*one(parent(A[1]))
      end
      return data(a)(vals[1:num]...)
   end
   if n < num
      A = vcat(A, zeros(Int, num - n))
   end
   return data(a)(A...)
end

function evaluate(a::UnivPoly{T}, vals::Vector{V}) where {T <: RingElement, V <: NCRingElem}
   return a(vals...)
end

function evaluate(a::UnivPoly{T}, vars::Vector{Int}, vals::Vector{V}) where {T <: RingElement, V <: RingElement}
   length(vars) != length(vals) && error("Numbers of variables and values do not match")
   vars2 = Vector{Int}(undef, 0)
   vals2 = Vector{mpoly_type(T)}(undef, 0)
   num = nvars(parent(data(a)))
   S = parent(a)
   n = nvars(S)
   for i = 1:length(vars)
      vars[i] > n && error("Unknown variable")
      if vars[i] <= num
         push!(vars2, vars[i])
         push!(vals2, data(S(vals[i])))
      end
   end
   return UnivPoly(evaluate(data(S(a)), vars2, vals2), S)
end

function evaluate(a::S, vars::Vector{S}, vals::Vector{V}) where {S <: UnivPoly{T}, V <: RingElement} where {T <: RingElement}
   varidx = Int[var_index(x) for x in vars]
   return evaluate(a, varidx, vals)
end

function (a::Union{MPolyRingElem, UniversalPolyRingElem})(;kwargs...)
   ss = symbols(parent(a))
   vars = Array{Int}(undef, length(kwargs))
   vals = Array{RingElement}(undef, length(kwargs))
   for (i, (var, val)) in enumerate(kwargs)
     vari = findfirst(isequal(var), ss)
     vari === nothing && error("Given polynomial has no variable $var")
     vars[i] = vari
     vals[i] = val
   end
   return evaluate(a, vars, vals)
end

########S,(a,b)=QQ[:a,:b]#######################################################################
#
#   GCD
#
###############################################################################

function gcd(a::UnivPoly{T}, b::UnivPoly{T}) where {T <: RingElement}
   a, b = univ_promote(a, b)
   return UnivPoly{T}(gcd(data(a), data(b)), parent(a))
end

function lcm(a::UnivPoly{T}, b::UnivPoly{T}) where {T <: RingElement}
   a, b = univ_promote(a, b)
   return UnivPoly{T}(lcm(data(a), data(b)), parent(a))
end

###############################################################################
#
#   Univariate polynomials
#
###############################################################################

function to_univariate(R::AbstractAlgebra.PolyRing{T}, p::UnivPoly{T}) where {T <: RingElement}
   return to_univariate(R, data(p))
end

to_univariate(p::UnivPoly) = to_univariate(data(p))

is_univariate(p::UnivPoly) = is_univariate(data(p))

is_univariate_with_data(p::UnivPoly) = is_univariate_with_data(data(p))

is_univariate(R::UniversalPolyRing) = is_univariate(mpoly_ring(R))

function coefficients_of_univariate(p::UnivPoly, check_univariate::Bool=true)
   return coefficients_of_univariate(data(p), check_univariate)
end

################################################################################
#
#  Change base ring
#
################################################################################

function _change_univ_poly_ring(R, Rx, cached::Bool)
   P = AbstractAlgebra.polynomial_ring_only(R, symbols(Rx); internal_ordering=internal_ordering(Rx), cached)
   S = universal_polynomial_ring(R; internal_ordering=internal_ordering(Rx), cached)
   S.mpoly_ring = P
   return S
end

function change_base_ring(R::Ring, p::UnivPoly{T}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(R, parent(p), cached)) where {T <: RingElement}
   base_ring(parent) != R && error("Base rings do not match.")
   return _map(R, p, parent)
end

function change_coefficient_ring(R::Ring, p::UnivPoly{T}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(R, parent(p), cached)) where {T <: RingElement}
  return change_base_ring(R, p, cached = cached, parent = parent)
end

################################################################################
#
#  Map
#
################################################################################

function map_coefficients(f::T, p::UnivPoly; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(parent(f(zero(base_ring(p)))), parent(p), cached)) where T
   return _map(f, p, parent)
end

function _map(g::T, p::UnivPoly, Rx) where T
   cvzip = zip(coefficients(p), exponent_vectors(p))
   M = MPolyBuildCtx(Rx)
   for (c, v) in cvzip
      push_term!(M, g(c), v)
   end

   return finish(M)
end

###############################################################################
#
#   MPolyBuildCtx
#
###############################################################################

function sort_terms!(p::UnivPoly{T}) where {T}
   p.p = sort_terms!(data(p))
   return p
end

function combine_like_terms!(p::UnivPoly{T}) where {T}
   p.p = combine_like_terms!(data(p))
   return p
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::AbstractAlgebra.UniversalPolyRing, _, _, _) = elem_type(S)

function RandomExtensions.make(S::AbstractAlgebra.UniversalPolyRing, term_range::AbstractUnitRange{Int},
                               exp_bound::AbstractUnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, exp_bound, vs[1])
   else
      Make(S, term_range, exp_bound, make(R, vs...))
   end
end

function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make4{
                 <:RingElement, <:AbstractAlgebra.UniversalPolyRing, <:AbstractUnitRange{Int}, <:AbstractUnitRange{Int}}})
   S, term_range, exp_bound, v = sp[][1:end]
   f = S()
   g = gens(S)
   R = base_ring(S)
   for i = 1:rand(rng, term_range)
      term = S(1)
      for j = 1:length(g)
         term *= g[j]^rand(rng, exp_bound)
      end
      term *= rand(rng, v)
      f += term
   end
   return f
end

function rand(rng::AbstractRNG, S::AbstractAlgebra.UniversalPolyRing,
              term_range::AbstractUnitRange{Int}, exp_bound::AbstractUnitRange{Int}, v...)
   rand(rng, make(S, term_range, exp_bound, v...))
end

function rand(S::AbstractAlgebra.UniversalPolyRing, term_range, exp_bound, v...)
   rand(Random.default_rng(), S, term_range, exp_bound, v...)
end

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::UniversalPolyRing{EuclideanRingResidueRingElem{BigInt}})
  return rand(R, 0:4, 0:10, -10:10)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::UnivPoly{T}) where {T <: RingElement}
  a.p = zero!(a.p)
  return a
end

function one!(a::UnivPoly{T}) where {T <: RingElement}
  a.p = one!(a.p)
  return a
end

function neg!(z::UnivPoly{T}, a::UnivPoly{T}) where {T <: RingElement}
  if parent(data(z)) == parent(data(a))
    z.p = neg!(z.p, a.p)
  else
    z.p = -a.p
  end
  return z
end

function fit!(a::UnivPoly, n::Int)
  fit!(data(a), n)
end

function add!(a::UnivPoly{T}, b::UnivPoly{T}, c::UnivPoly{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = add!(data(a), data(b), data(c))
  else
    a.p = data(b + c)
  end
  return a
end

function add!(a::UnivPoly{T}, b::UnivPoly{T}, c::RingElement) where {T <: RingElement}
  if parent(data(a)) == parent(data(b))
    a.p = add!(data(a), data(b), c)
  else
    a.p = data(b + c)
  end
  return a
end

add!(a::UnivPoly{T}, b::RingElement, c::UnivPoly{T}) where {T <: RingElement} = add!(a, c, b)

function sub!(a::UnivPoly{T}, b::UnivPoly{T}, c::UnivPoly{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = sub!(data(a), data(b), data(c))
  else
    a.p = data(b - c)
  end
  return a
end

function sub!(a::UnivPoly{T}, b::UnivPoly{T}, c::RingElement) where {T <: RingElement}
  if parent(data(a)) == parent(data(b))
    a.p = sub!(data(a), data(b), c)
  else
    a.p = data(b - c)
  end
  return a
end

function sub!(a::UnivPoly{T}, b::RingElement, c::UnivPoly{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(c))
    a.p = sub!(data(a), b, data(c))
  else
    a.p = data(b - c)
  end
  return a
end

function mul!(a::UnivPoly{T}, b::UnivPoly{T}, c::UnivPoly{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = mul!(data(a), data(b), data(c))
  else
    a.p = data(b * c)
  end
  return a
end

function mul!(a::UnivPoly{T}, b::UnivPoly{T}, c::RingElement) where {T <: RingElement}
  if parent(data(a)) == parent(data(b))
    a.p = mul!(data(a), data(b), c)
  else
    a.p = data(b * c)
  end
  return a
end

mul!(a::UnivPoly{T}, b::RingElement, c::UnivPoly{T}) where {T <: RingElement} = mul!(a, c, b)

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{UnivPoly{T}}, ::Type{UnivPoly{T}}) where {T <: RingElement} = UnivPoly{T}

function promote_rule(::Type{UnivPoly{T}}, ::Type{V}) where {T <: RingElement, V <: RingElement}
   promote_rule(T, V) == T ? UnivPoly{T} : Union{}
end

###############################################################################
#
#   Parent object overload
#
###############################################################################

function upgrade(S::UniversalPolyRing{T}, pp::MPolyRingElem{T}) where {T}
   n = nvars(S) - nvars(parent(pp))
   ctx = MPolyBuildCtx(mpoly_ring(S))
   v0 = zeros(Int, n)
   for (c, v) in zip(coefficients(pp), exponent_vectors(pp))
      push_term!(ctx, c, vcat(v, v0))
   end
   return finish(ctx)
end

function (a::UniversalPolyRing{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::UniversalPolyRing{T})() where {T <: RingElement}
   return UnivPoly{T}(mpoly_ring(a)(), a)
end

function (a::UniversalPolyRing{T})(b::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
   return UnivPoly{T}(mpoly_ring(a)(b), a)
end

function (a::UniversalPolyRing{T})(b::T) where {T <: RingElem}
   return UnivPoly{T}(mpoly_ring(a)(b), a)
end

function (S::UniversalPolyRing{T})(p::UnivPoly{T}) where {T <: RingElement}
   parent(p) !== S && error("Unable to coerce")
   n = nvars(S) - nvars(parent(data(p)))
   if n != 0
      p = UnivPoly{T}(upgrade(S, data(p)), S)
   end
   return p
end

function (a::UniversalPolyRing{T})(b::Vector{T}, m::Vector{Vector{Int}}) where {T <: RingElement}
   if length(m) != 0
      len = length(m[1])
      num = nvars(mpoly_ring(a))
      if len != num
         for i = 1:length(m)
            m[i] = vcat(m[i], zeros(Int, num - len))
         end
      end
   end
   return UnivPoly{T}(mpoly_ring(a)(b, m), a)
end
