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

base_ring_type(::Type{UniversalRing{T}}) where T = parent_type(T)
base_ring(S::UniversalRing) = S.base_ring::base_ring_type(S)

coefficient_ring_type(::Type{<:UniversalRing{<:MPolyRingElem{T}}}) where T = parent_type(T)
coefficient_ring(S::UniversalRing{<:MPolyRingElem}) = coefficient_ring(base_ring(S))::coefficient_ring_type(S)

is_domain_type(::Type{UniversalRingElem{T}}) where T = is_domain_type(T)

is_exact_type(::Type{UniversalRingElem{T}}) where T = is_exact_type(T)

parent(p::UniversalRingElem) = p.parent::parent_type(p)

elem_type(::Type{UniversalRing{T}}) where T = UniversalRingElem{T}

parent_type(::Type{UniversalRingElem{T}}) where T = UniversalRing{T}

number_of_variables(S::UniversalRing) = number_of_variables(base_ring(S))

number_of_generators(S::UniversalRing) = number_of_generators(base_ring(S))

symbols(S::UniversalRing) = symbols(base_ring(S))

internal_ordering(p::UniversalPolyRing) = internal_ordering(base_ring(p))

data(p::UniversalRingElem) = p.p

###############################################################################
#
#   Manipulating terms and monomials
#
###############################################################################

exponent_vector(p::UniversalRingElem{<:MPolyRingElem}, i::Int) = exponent_vector(data(p), i)

function exponent(p::UniversalRingElem{<:MPolyRingElem}, i::Int, j::Int)
   if j > nvars(parent(data(p))) && j <= nvars(parent(p))
      return 0
   end
   return exponent(data(p), i, j)
end

function set_exponent_vector!(p::UniversalRingElem{<:MPolyRingElem}, i::Int, exps::Vector{Int})
   S = parent(p)
   len = length(exps)
   upgrade!(p)
   if len < nvars(S)
      exps = vcat(exps, zeros(Int, nvars(S) - len))
   end
   p.p = set_exponent_vector!(data(p), i, exps)
   return p
end

function coeff(p::UniversalRingElem{<:MPolyRingElem}, exps::Vector{Int})
   S = parent(p)
   len = length(exps)
   n = nvars(parent(data(p)))
   if len > n
      if !iszero(exps[n + 1:len])
         return coefficient_ring(S)()
      end
      return coeff(data(p), exps[1:n])
   end
   n = nvars(parent(data(p)))
   if len < n
      exps = vcat(exps, zeros(Int, n - len))
   end
   return coeff(data(p), exps)
end

function setcoeff!(p::UniversalRingElem{<:MPolyRingElem}, exps::Vector{Int}, c::RingElement)
   c = coefficient_ring(data(p))(c)
   S = parent(p)
   len = length(exps)
   upgrade!(p)
   if len < nvars(S)
      exps = vcat(exps, zeros(Int, nvars(S) - len))
   end
   p.p = setcoeff!(data(p), exps, c)
   return p
end

function setcoeff!(a::UniversalRingElem{<:MPolyRingElem}, i::Int, c::RingElement)
   setcoeff!(data(a), i, c)
   return a
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(R::UniversalRing{T}) where T = UniversalRingElem{T}(zero(base_ring(R)), R)

one(R::UniversalRing{T}) where T = UniversalRingElem{T}(one(base_ring(R)), R)

iszero(p::UniversalRingElem) = iszero(data(p))

isone(p::UniversalRingElem) = isone(data(p))

is_unit(p::UniversalRingElem) = is_unit(data(p))

is_zero_divisor(p::UniversalRingElem) = is_zero_divisor(data(p))

function is_zero_divisor_with_annihilator(p::UniversalRingElem{T}) where T
   f, b = is_zero_divisor_with_annihilator(data(p))
   return f, UniversalRingElem{T}(b, parent(p))
end

is_gen(p::UniversalRingElem) = is_gen(data(p))

is_homogeneous(p::UniversalRingElem{<:MPolyRingElem}) = is_homogeneous(data(p))

is_monomial(p::UniversalRingElem{<:MPolyRingElem}) = is_monomial(data(p))

is_constant(p::UniversalRingElem{<:MPolyRingElem}) = is_constant(data(p))

is_term(p::UniversalRingElem{<:MPolyRingElem}) = is_term(data(p))

coeff(p::UniversalRingElem{<:MPolyRingElem}, i::Int) = coeff(data(p), i)

function coeff(p::UniversalRingElem{<:MPolyRingElem{T}}, m::UniversalRingElem{<:MPolyRingElem{T}}) where T
   p, m = univ_promote(p, m)
   return coeff(data(p), data(m))
end

function monomial(p::UniversalRingElem{<:MPolyRingElem{T}}, i::Int) where T
   S = parent(p)
   m = monomial(data(p), i)
   return UniversalRingElem(m, S)
end

function monomial!(m::UniversalRingElem{<:MPolyRingElem{T}}, p::UniversalRingElem{<:MPolyRingElem{T}}, i::Int) where T
   parent(m) != parent(p) && error("Incompatible monomial")
   if parent(data(m)) != parent(data(p))
      m.p = parent(data(p))()
   end
   m.p = monomial!(data(m), data(p), i)
   return m
end

function term(p::UniversalRingElem{<:MPolyRingElem{T}}, i::Int) where T
   S = parent(p)
   t = term(data(p), i)
   return UniversalRingElem(t, S)
end

leading_coefficient(p::UniversalRingElem{<:MPolyRingElem}) = leading_coefficient(data(p))

trailing_coefficient(p::UniversalRingElem{<:MPolyRingElem}) = trailing_coefficient(data(p))

function tail(p::UniversalRingElem{<:MPolyRingElem{T}}) where T
   S = parent(p)
   return UniversalRingElem(tail(data(p)), S)
end

constant_coefficient(p::UniversalRingElem{<:MPolyRingElem}) = constant_coefficient(data(p))

function leading_monomial(p::UniversalRingElem{<:MPolyRingElem{T}}) where T
   S = parent(p)
   return UniversalRingElem(leading_monomial(data(p)), S)
end

function leading_term(p::UniversalRingElem{<:MPolyRingElem{T}}) where T
   S = parent(p)
   return UniversalRingElem(leading_term(data(p)), S)
end

max_fields(p::UniversalRingElem{<:MPolyRingElem}) = max_fields(data(p))

function degree(p::UniversalRingElem{<:MPolyRingElem}, i::Int)
   if i <= nvars(parent(p)) && i > nvars(parent(data(p)))
      return 0
   end
   return degree(data(p), i)
end

function degree(f::UniversalRingElem{<:MPolyRingElem{T}}, x::UniversalRingElem{<:MPolyRingElem{T}}) where T
   check_parent(f, x)
   return degree(f, var_index(x))
end

function degrees(p::UniversalRingElem{<:MPolyRingElem})
   v = degrees(data(p))
   len = length(v)
   num = nvars(parent(p))
   if len < num
      v = vcat(v, zeros(Int, num - len))
   end
   return v
end

total_degree(p::UniversalRingElem{<:MPolyRingElem}) = total_degree(data(p))

length(p::UniversalRingElem{<:MPolyRingElem}) = length(data(p))

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

function gen(S::UniversalRing{T}, i::Int) where T
   @boundscheck 1 <= i <= nvars(S) || throw(ArgumentError("generator index out of range"))
   return UniversalRingElem{T}(gen(base_ring(S), i), S)
end

function gens(S::UniversalRing{T}) where T
   return [UniversalRingElem{T}(g, S) for g in gens(base_ring(S))]
end

# HACK: we abuse the @varnames_interface macro to teach gens for UniversalRing
# some superpowers
function _univ_gens(S::UniversalRing{T}, vars::Vector{Symbol}) where T
   idx = _ensure_variables(S, vars)
   # TRICK: @varnames_interface expects two return values, but we only care
   # for the second; so just return literally nothing for the first
   return nothing, [UniversalRingElem{T}(gen(base_ring(S), i), S) for i in idx]
end

AbstractAlgebra.@varnames_interface _univ_gens(R::UniversalRing{T}, s) where T

function gens(S::UniversalRing{T}, a::AbstractAlgebra.VarNames, as::AbstractAlgebra.VarNames...) where T
   res = _univ_gens(S, a, as...)
   length(res) == 2 && return res[2] # special case for improved backwards compatibility
   return res[2:end]
end

var_index(x::UniversalRingElem{<:MPolyRingElem}) = var_index(data(x))

var_indices(p::UnivPoly) = var_indices(data(p))

function vars(p::UniversalRingElem{<:MPolyRingElem{T}}) where {T <: RingElement}
   V = vars(data(p))
   S = parent(p)
   return [UniversalRingElem(v, S) for v in V]
end

canonical_unit(p::UniversalRingElem) = canonical_unit(data(p))

characteristic(R::UniversalRing) = characteristic(base_ring(R))
is_known(::typeof(characteristic), R::UniversalRing) = is_known(characteristic, base_ring(R))

function Base.hash(p::UniversalRingElem, h::UInt)  # TODO
   b = 0xcf418d4529109236%UInt
   return xor(b, hash(data(p), h))
end

function deepcopy_internal(p::UniversalRingElem{T}, dict::IdDict) where T
   return UniversalRingElem{T}(deepcopy_internal(data(p), dict), parent(p))
end

Base.copy(f::UniversalRingElem{T}) where T = UniversalRingElem{T}(copy(data(f)), parent(f))

###############################################################################
#
#   Multivariate coefficients
#
###############################################################################

function coeff(p::UniversalRingElem{<:MPolyRingElem{T}}, vars::Vector{Int}, exps::Vector{Int}) where {T}
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
   return UniversalRingElem(coeff(data(p), vars2, exps2), S)
end

function coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where {U, T <: UniversalRingElem{<:MPolyRingElem{U}}}
   varidx = [var_index(x) for x in vars]
   return coeff(a, varidx, exps)
end

###############################################################################
#
#   String I/O
#
###############################################################################

universal_ring_name(R::UniversalRing) = "ring"

universal_ring_name(R::UniversalRing{<:MPolyRingElem}) = "polynomial ring"

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

function -(p::UniversalRingElem{T}) where T
   return UniversalRingElem{T}(-data(p), parent(p))
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

function Base.sqrt(p::UniversalRingElem{T}; check::Bool=true) where T
   S = parent(p)
   s = sqrt(data(p); check=check)
   return UniversalRingElem{T}(s, S)
end

is_square(p::UniversalRingElem) = is_square(data(p))

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

for op in (:+, :-, :*)
  @eval begin
    function $op(p::UniversalRingElem{T}, n::JuliaRingElement) where {T}
       S = parent(p)
       return UniversalRingElem{T}($op(data(p),n), S)
    end

    function $op(p::UniversalRingElem{T}, n::T) where {T <: RingElem}
       S = parent(p)
       return UniversalRingElem{T}($op(data(p),n), S)
    end

    function $op(n::JuliaRingElement, p::UniversalRingElem{T}) where {T}
       S = parent(p)
       return UniversalRingElem{T}($op(n,data(p)), S)
    end

    function $op(n::T, p::UniversalRingElem{T}) where {T <: RingElem}
       S = parent(p)
       return UniversalRingElem{T}($op(n,data(p)), S)
    end
  end
end

function divexact(p::UniversalRingElem{T}, n::JuliaRingElement; check::Bool=true) where T
   S = parent(p)
   return UniversalRingElem{T}(divexact(data(p), n; check=check), S)
end

function divexact(p::UniversalRingElem{T}, n::T; check::Bool=true) where {T <: RingElem}
   S = parent(p)
   return UniversalRingElem{T}(divexact(data(p), n; check=check), S)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::UniversalRingElem{T}, b::UniversalRingElem{T}) where {T}
   check_parent(a, b)
   upgrade!(a)
   upgrade!(b)
   return data(a) == data(b)
end

function ==(a::UniversalRingElem{<:MPolyRingElem{T}}, b::UniversalRingElem{<:MPolyRingElem{T}}) where {T}
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

function isless(a::UniversalRingElem{<:MPolyRingElem{T}}, b::UniversalRingElem{<:MPolyRingElem{T}}) where {T}
   check_parent(a, b)
   return isless(data(upgrade!(a)), data(upgrade!(b)))
end

###############################################################################
#
#   Ad hoc comparison functions
#
###############################################################################

==(p::UniversalRingElem, n::RingElement) = data(p) == n

==(n::RingElement, p::UniversalRingElem) = data(p) == n

###############################################################################
#
#   Powering
#
###############################################################################

function ^(p::UniversalRingElem{T}, b::Int) where T
   S = parent(p)
   return UniversalRingElem{T}(data(p)^b, S)
end

###############################################################################
#
#   Inflation/deflation
#
###############################################################################

deflation(p::UniversalRingElem{<:MPolyRingElem{T}}) where T = deflation(data(p))

function deflate(p::UniversalRingElem{<:MPolyRingElem{T}}, shift::Vector{Int}, defl::Vector{Int}) where T
   S = parent(p)
   num = nvars(S)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   upgrade!(p)
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return UniversalRingElem(deflate(data(p), shift, defl), S)
end

function deflate(p::UniversalRingElem{<:MPolyRingElem{T}}, defl::Vector{Int}) where T
   return deflate(p, zeros(Int, length(defl)), defl)
end

function inflate(p::UniversalRingElem{<:MPolyRingElem{T}}, shift::Vector{Int}, defl::Vector{Int}) where T
   S = parent(p)
   num = nvars(S)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   upgrade!(p)
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return UniversalRingElem(inflate(data(p), shift, defl), S)
end

function inflate(p::UniversalRingElem{<:MPolyRingElem{T}}, defl::Vector{Int}) where T
   return inflate(p, zeros(Int, length(defl)), defl)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::UniversalRingElem{T}, b::UniversalRingElem{T}; check::Bool=true) where T
   a, b = univ_promote(a, b)
   return UniversalRingElem{T}(divexact(data(a), data(b); check=check), parent(a))
end

function divides(a::UniversalRingElem{T}, b::UniversalRingElem{T}) where T
   a, b = univ_promote(a, b)
   flag, q = divides(data(a), data(b))
   return flag, UniversalRingElem{T}(q, parent(a))
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function Base.div(a::UniversalRingElem{T}, b::UniversalRingElem{T}) where T
   a, b = univ_promote(a, b)
   return UniversalRingElem{T}(div(data(a), data(b)), parent(a))
end

function Base.divrem(a::UniversalRingElem{T}, b::UniversalRingElem{T}) where T
   a, b = univ_promote(a, b)
   q, r = divrem(data(a), data(b))
   return UniversalRingElem{T}(q, parent(a)), UniversalRingElem{T}(r, parent(a))
end

###############################################################################
#
#   Derivative
#
###############################################################################

function derivative(p::UniversalRingElem{<:MPolyRingElem{T}}, j::Int) where {T}
   j > nvars(parent(p)) && error("No such variable")
   if j > nvars(parent(data(p)))
      return zero(parent(p))
   end
   return UniversalRingElem(derivative(data(p), j), parent(p))
end

function derivative(p::UniversalRingElem{<:MPolyRingElem{T}}, x::UniversalRingElem{<:MPolyRingElem{T}}) where {T}
   return derivative(p, var_index(x))
end

###############################################################################
#
#   Remove and valuation
#
###############################################################################

function remove(z::UniversalRingElem{<:MPolyRingElem{T}}, p::UniversalRingElem{<:MPolyRingElem{T}}) where {T}
   z, p = univ_promote(z, p)
   val, q = remove(data(z), data(p))
   return val, UniversalRingElem(q, parent(z))
end

function valuation(z::UniversalRingElem{<:MPolyRingElem{T}}, p::UniversalRingElem{<:MPolyRingElem{T}}) where {T}
  v, _ = remove(z, p)
  return v
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(a::UnivPoly, A::Vector{<:NCRingElement})
   isempty(A) && error("Evaluating at an empty list of values is not allowed")
   a2 = data(a)
   varidx = var_indices(a2)
   isempty(varidx) && return constant_coefficient(a2) # TODO: this is weird
   vals = [zero(parent(A[1])) for _ in 1:nvars(parent(a2))]
   n = length(A)
   for i in varidx
      i <= n || error("Number of variables does not match number of values")
      vals[i] = A[i]
   end
   return evaluate(a2, vals)
end

function (a::UnivPoly)(val::T, vals::T...) where T <: NCRingElement
   return evaluate(a, [val, vals...])
end

function (a::UnivPoly)()
   return evaluate(a, Int[])
end

function evaluate(a::UniversalRingElem{<:MPolyRingElem{T}}, vars::Vector{Int}, vals::Vector{V}) where {T <: RingElement, V <: RingElement}
   length(vars) != length(vals) && error("Numbers of variables and values do not match")
   vars2 = Vector{Int}(undef, 0)
   vals2 = Vector{mpoly_type(T)}(undef, 0)
   num = nvars(parent(data(a)))
   S = parent(a)
   upgrade!(a)
   a2 = data(a)
   for i = 1:length(vars)
      if vars[i] <= num
         push!(vars2, vars[i])
         push!(vals2, data(S(vals[i])))
      end
   end
   return UniversalRingElem(evaluate(a2, vars2, vals2), S)
end

function evaluate(a::S, vars::Vector{S}, vals::Vector{V}) where {S <: UniversalRingElem{<:MPolyRingElem{T}}, V <: RingElement} where {T <: RingElement}
   varidx = Int[var_index(x) for x in vars]
   return evaluate(a, varidx, vals)
end

########S,(a,b)=QQ[:a,:b]#######################################################################
#
#   GCD
#
###############################################################################

function gcd(a::UniversalRingElem{T}, b::UniversalRingElem{T}) where {T <: RingElement}
   a, b = univ_promote(a, b)
   return UniversalRingElem{T}(gcd(data(a), data(b)), parent(a))
end

function lcm(a::UniversalRingElem{T}, b::UniversalRingElem{T}) where {T <: RingElement}
   a, b = univ_promote(a, b)
   return UniversalRingElem{T}(lcm(data(a), data(b)), parent(a))
end

###############################################################################
#
#   Univariate polynomials
#
###############################################################################

function to_univariate(R::AbstractAlgebra.PolyRing{T}, p::UniversalRingElem{<:MPolyRingElem{T}}) where {T <: RingElement}
   return to_univariate(R, data(p))
end

to_univariate(p::UniversalRingElem{<:MPolyRingElem}) = to_univariate(data(p))

is_univariate(p::UniversalRingElem{<:MPolyRingElem}) = is_univariate(data(p))

is_univariate_with_data(p::UniversalRingElem{<:MPolyRingElem}) = is_univariate_with_data(data(p))

is_univariate(R::UniversalRing{<:MPolyRingElem}) = is_univariate(base_ring(R))

function coefficients_of_univariate(p::UniversalRingElem{<:MPolyRingElem}, check_univariate::Bool=true)
   return coefficients_of_univariate(data(p), check_univariate)
end

################################################################################
#
#  Change base ring
#
################################################################################

_change_univ_poly_ring(R, Rx, cached::Bool) = universal_polynomial_ring(R, symbols(Rx); internal_ordering=internal_ordering(Rx), cached)[1]

function change_base_ring(R::Ring, p::UniversalRingElem{<:MPolyRingElem{T}}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(R, parent(p), cached)) where {T <: RingElement}
   upgrade!(p)
   return UniversalRingElem(change_base_ring(R, data(p); parent = base_ring(parent)), parent)
end

function change_coefficient_ring(R::Ring, p::UniversalRingElem{<:MPolyRingElem{T}}; cached::Bool=true, parent::UniversalRing{<:MPolyRingElem} = _change_univ_poly_ring(R, parent(p), cached)) where {T <: RingElement}
  return change_base_ring(R, p, cached = cached, parent = parent)
end

################################################################################
#
#  Map
#
################################################################################

function map_coefficients(f::Any, p::UniversalRingElem{<:MPolyRingElem{T}}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(parent(f(zero(coefficient_ring(p)))), parent(p), cached)) where T <: RingElement
   upgrade!(p)
   return UniversalRingElem(map_coefficients(f, data(p); parent = base_ring(parent)), parent)
end

###############################################################################
#
#   MPolyBuildCtx
#
###############################################################################

function sort_terms!(p::UniversalRingElem{<:MPolyRingElem})
   p.p = sort_terms!(data(p))
   return p
end

function combine_like_terms!(p::UniversalRingElem{<:MPolyRingElem})
   p.p = combine_like_terms!(data(p))
   return p
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
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::UniversalRing{MPoly{EuclideanRingResidueRingElem{BigInt}}})
  return rand(R, 0:4, 0:10, -10:10)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::UniversalRingElem{T}) where {T <: RingElement}
  a.p = zero!(a.p)
  return a
end

function one!(a::UniversalRingElem{T}) where {T <: RingElement}
  a.p = one!(a.p)
  return a
end

function neg!(z::UniversalRingElem{T}, a::UniversalRingElem{T}) where {T <: RingElement}
  if parent(data(z)) == parent(data(a))
    z.p = neg!(z.p, a.p)
  else
    z.p = -a.p
  end
  return z
end

function fit!(a::UniversalRingElem{<:MPolyRingElem}, n::Int)
  fit!(data(a), n)
end

function add!(a::UniversalRingElem{T}, b::UniversalRingElem{T}, c::UniversalRingElem{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = add!(data(a), data(b), data(c))
  else
    a.p = data(b + c)
  end
  return a
end

function add!(a::UniversalRingElem{T}, b::UniversalRingElem{T}, c::RingElement) where {T <: RingElement}
  if parent(data(a)) == parent(data(b))
    a.p = add!(data(a), data(b), c)
  else
    a.p = data(b + c)
  end
  return a
end

add!(a::UniversalRingElem{T}, b::RingElement, c::UniversalRingElem{T}) where {T <: RingElement} = add!(a, c, b)

function sub!(a::UniversalRingElem{T}, b::UniversalRingElem{T}, c::UniversalRingElem{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = sub!(data(a), data(b), data(c))
  else
    a.p = data(b - c)
  end
  return a
end

function sub!(a::UniversalRingElem{T}, b::UniversalRingElem{T}, c::RingElement) where {T <: RingElement}
  if parent(data(a)) == parent(data(b))
    a.p = sub!(data(a), data(b), c)
  else
    a.p = data(b - c)
  end
  return a
end

function sub!(a::UniversalRingElem{T}, b::RingElement, c::UniversalRingElem{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(c))
    a.p = sub!(data(a), b, data(c))
  else
    a.p = data(b - c)
  end
  return a
end

function mul!(a::UniversalRingElem{T}, b::UniversalRingElem{T}, c::UniversalRingElem{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = mul!(data(a), data(b), data(c))
  else
    a.p = data(b * c)
  end
  return a
end

function mul!(a::UniversalRingElem{T}, b::UniversalRingElem{T}, c::RingElement) where {T <: RingElement}
  if parent(data(a)) == parent(data(b))
    a.p = mul!(data(a), data(b), c)
  else
    a.p = data(b * c)
  end
  return a
end

mul!(a::UniversalRingElem{T}, b::RingElement, c::UniversalRingElem{T}) where {T <: RingElement} = mul!(a, c, b)

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{UniversalRingElem{T}}, ::Type{UniversalRingElem{T}}) where {T <: RingElement} = UniversalRingElem{T}

function promote_rule(::Type{UniversalRingElem{T}}, ::Type{V}) where {T <: RingElement, V <: RingElement}
   promote_rule(T, V) == T ? UniversalRingElem{T} : Union{}
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

function (a::UniversalRing)(b::RingElement)
   return a(base_ring(a)(b))
end

function (a::UniversalRing{T})() where {T <: RingElement}
   return UniversalRingElem{T}(base_ring(a)(), a)
end

function (a::UniversalRing{T})(b::T) where {T <: RingElem}
   return UniversalRingElem{T}(base_ring(a)(b), a)
end

function (S::UniversalRing{T})(p::UniversalRingElem{T}) where {T <: RingElement}
   parent(p) !== S && error("Unable to coerce")
   return upgrade!(p)
end

function (a::UniversalRing{<:MPolyRingElem{T}})(b::Vector{T}, m::Vector{Vector{Int}}) where {T <: RingElement}
   if length(m) != 0
      len = length(m[1])
      num = nvars(base_ring(a))
      if len != num
         for i = 1:length(m)
            m[i] = vcat(m[i], zeros(Int, num - len))
         end
      end
   end
   return UniversalRingElem(base_ring(a)(b, m), a)
end
