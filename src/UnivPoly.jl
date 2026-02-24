###############################################################################
#
#   UnivPoly.jl : Universal polynomial ring (variables can be added)
#
###############################################################################

internal_ordering(p::UniversalRing{<:MPolyRingElem}) = internal_ordering(base_ring(p))

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

is_homogeneous(p::UniversalRingElem{<:MPolyRingElem}) = is_homogeneous(data(p))

is_monomial(p::UniversalRingElem{<:MPolyRingElem}) = is_monomial(data(p))

is_constant(p::UniversalRingElem{<:MPolyRingElem}) = is_constant(data(p))

is_term(p::UniversalRingElem{<:MPolyRingElem}) = is_term(data(p))

coeff(p::UniversalRingElem{<:MPolyRingElem}, i::Int) = coeff(data(p), i)

function coeff(p::T, m::T) where {T <: UniversalRingElem{<:MPolyRingElem}}
   p, m = univ_promote(p, m)
   return coeff(data(p), data(m))
end

function content(a::UniversalRingElem{<:MPolyRingElem})
   return content(data(a))
end

function monomial(p::UniversalRingElem{<:MPolyRingElem}, i::Int)
   S = parent(p)
   m = monomial(data(p), i)
   return UniversalRingElem(m, S)
end

function monomial!(m::T, p::T, i::Int) where {T <: UniversalRingElem{<:MPolyRingElem}}
   parent(m) != parent(p) && error("Incompatible monomial")
   if parent(data(m)) != parent(data(p))
      m.p = parent(data(p))()
   end
   m.p = monomial!(data(m), data(p), i)
   return m
end

function term(p::UniversalRingElem{<:MPolyRingElem}, i::Int)
   S = parent(p)
   t = term(data(p), i)
   return UniversalRingElem(t, S)
end

leading_coefficient(p::UniversalRingElem{<:MPolyRingElem}) = leading_coefficient(data(p))

trailing_coefficient(p::UniversalRingElem{<:MPolyRingElem}) = trailing_coefficient(data(p))

function tail(p::UniversalRingElem{<:MPolyRingElem})
   S = parent(p)
   return UniversalRingElem(tail(data(p)), S)
end

constant_coefficient(p::UniversalRingElem{<:MPolyRingElem}) = constant_coefficient(data(p))

function leading_monomial(p::UniversalRingElem{<:MPolyRingElem})
   S = parent(p)
   return UniversalRingElem(leading_monomial(data(p)), S)
end

function leading_term(p::UniversalRingElem{<:MPolyRingElem})
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

function degree(f::T, x::T) where {T <: UniversalRingElem{<:MPolyRingElem}}
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

var_index(x::UniversalRingElem{<:MPolyRingElem}) = var_index(data(x))

var_indices(p::UniversalRingElem{<:MPolyRingElem}) = var_indices(data(p))

function vars(p::UniversalRingElem{<:MPolyRingElem})
   V = vars(data(p))
   S = parent(p)
   return [UniversalRingElem(v, S) for v in V]
end

function Base.hash(p::UniversalRingElem{<:MPolyRingElem}, h::UInt)
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

###############################################################################
#
#   Multivariate coefficients
#
###############################################################################

function coeff(p::UniversalRingElem{<:MPolyRingElem}, vars::Vector{Int}, exps::Vector{Int})
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

function coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where {T <: UniversalRingElem{<:MPolyRingElem}}
   varidx = [var_index(x) for x in vars]
   return coeff(a, varidx, exps)
end

###############################################################################
#
#   String I/O
#
###############################################################################

universal_ring_name(R::UniversalRing{<:MPolyRingElem}) = "polynomial ring"

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

function Base.eltype(::Type{UnivPolyCoeffs{UniversalRingElem{<:MPolyRingElem, T}}}) where T <: RingElement
   return T
end

function Base.eltype(::Type{UnivPolyExponentVectors{T}}) where T <: UniversalRingElem{<:MPolyRingElem}
   return Vector{Int}
end

function Base.eltype(::Type{UnivPolyMonomials{T}}) where T <: UniversalRingElem{<:MPolyRingElem}
   return T
end

function Base.eltype(::Type{UnivPolyTerms{T}}) where T <: UniversalRingElem{<:MPolyRingElem}
   return T
end

function coefficients(a::UniversalRingElem{<:MPolyRingElem})
   return UnivPolyCoeffs(a)
end

function exponent_vectors(a::UniversalRingElem{<:MPolyRingElem})
   return UnivPolyExponentVectors(a)
end

function monomials(a::UniversalRingElem{<:MPolyRingElem})
   return UnivPolyMonomials(a)
end

function terms(a::UniversalRingElem{<:MPolyRingElem})
   return UnivPolyTerms(a)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::T, b::T) where {T <: UniversalRingElem{<:MPolyRingElem}}
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

function isless(a::T, b::T) where {T <: UniversalRingElem{<:MPolyRingElem}}
   check_parent(a, b)
   return isless(data(upgrade!(a)), data(upgrade!(b)))
end

###############################################################################
#
#   Inflation/deflation
#
###############################################################################

deflation(p::UniversalRingElem{<:MPolyRingElem}) = deflation(data(p))

function deflate(p::T, shift::Vector{Int}, defl::Vector{Int}) where {T <: UniversalRingElem{<:MPolyRingElem}}
   S = parent(p)
   num = nvars(S)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   upgrade!(p)
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return T(deflate(data(p), shift, defl), S)
end

function deflate(p::UniversalRingElem{<:MPolyRingElem}, defl::Vector{Int})
   return deflate(p, zeros(Int, length(defl)), defl)
end

function inflate(p::T, shift::Vector{Int}, defl::Vector{Int}) where {T <: UniversalRingElem{<:MPolyRingElem}}
   S = parent(p)
   num = nvars(S)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   upgrade!(p)
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return T(inflate(data(p), shift, defl), S)
end

function inflate(p::UniversalRingElem{<:MPolyRingElem}, defl::Vector{Int})
   return inflate(p, zeros(Int, length(defl)), defl)
end

###############################################################################
#
#   Derivative
#
###############################################################################

function derivative(p::T, j::Int) where {T <: UniversalRingElem{<:MPolyRingElem}}
   j > nvars(parent(p)) && error("No such variable")
   if j > nvars(parent(data(p)))
      return zero(parent(p))
   end
   return T(derivative(data(p), j), parent(p))
end

function derivative(p::T, x::T) where {T <: UniversalRingElem{<:MPolyRingElem}}
   return derivative(p, var_index(x))
end

###############################################################################
#
#   Remove and valuation
#
###############################################################################

function remove(z::T, p::T) where {T <: UniversalRingElem{<:MPolyRingElem}}
   z, p = univ_promote(z, p)
   val, q = remove(data(z), data(p))
   return val, UniversalRingElem(q, parent(z))
end

function valuation(z::T, p::T) where {T <: UniversalRingElem{<:MPolyRingElem}}
  v, _ = remove(z, p)
  return v
end

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

@doc raw"""
    universal_poly_type(::Type{T}) where T<:RingElement
    universal_poly_type(::T) where T<:RingElement
    universal_poly_type(::Type{S}) where S<:Ring
    universal_poly_type(::S) where S<:Ring

The type of universal polynomials with coefficients of type `T` respectively `elem_type(S)`.

See also [`universal_poly_ring_type`](@ref), [`mpoly_type`](@ref) and [`mpoly_ring_type`](@ref).
"""
universal_poly_type(::Type{T}) where T<:RingElement = UniversalRingElem{mpoly_type(T), T}
universal_poly_type(::Type{S}) where S<:Ring = universal_poly_type(elem_type(S))
universal_poly_type(x) = universal_poly_type(typeof(x)) # to stop this method from eternally recursing on itself, we better add ...
universal_poly_type(::Type{T}) where T = throw(ArgumentError("Type `$T` must be subtype of `RingElement`."))
universal_poly_type(T::Type{Union{}}) = throw(MethodError(universal_poly_type, (T,)))

@doc raw"""
    universal_poly_ring_type(::Type{T}) where T<:RingElement
    universal_poly_ring_type(::T) where T<:RingElement
    universal_poly_ring_type(::Type{S}) where S<:Ring
    universal_poly_ring_type(::S) where S<:Ring

The type of universal polynomial rings with coefficients of type `T`
respectively `elem_type(S)`. Implemented via [`universal_poly_type`](@ref).

See also [`mpoly_type`](@ref) and [`mpoly_ring_type`](@ref).
"""
universal_poly_ring_type(x) = parent_type(universal_poly_type(x))

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(a::UniversalRingElem{<:MPolyRingElem}, A::Vector{<:NCRingElement})
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

function (a::UniversalRingElem{<:MPolyRingElem})(val::T, vals::T...) where T <: NCRingElement
   return evaluate(a, [val, vals...])
end

function evaluate(a::UniversalRingElem{T, U}, vars::Vector{Int}, vals::Vector{<:RingElement}) where {T <: MPolyRingElem, U}
   length(vars) != length(vals) && error("Numbers of variables and values do not match")
   vars2 = Vector{Int}(undef, 0)
   vals2 = Vector{T}(undef, 0)
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
   return UniversalRingElem{T, U}(evaluate(a2, vars2, vals2), S)
end

function evaluate(a::T, vars::Vector{T}, vals::Vector{<:RingElement}) where {T <: UniversalRingElem{<:MPolyRingElem}}
   varidx = Int[var_index(x) for x in vars]
   return evaluate(a, varidx, vals)
end

function (a::UniversalRingElem{<:MPolyRingElem})(;kwargs...)
   ss = symbols(parent(a))
   vars = Int[]
   vals = RingElement[]
   for (var, val) in kwargs
     vari = findfirst(isequal(var), ss)
     vari === nothing && continue
     push!(vars, vari)
     push!(vals, val)
   end
   return evaluate(a, vars, vals)
end

###############################################################################
#
#   Univariate polynomials
#
###############################################################################

function to_univariate(R::AbstractAlgebra.PolyRing{T}, p::UniversalRingElem{<:MPolyRingElem, T}) where T
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

function change_base_ring(R::Ring, p::UniversalRingElem{<:MPolyRingElem, T}; cached::Bool=true, parent::UniversalRing{<:MPolyRingElem} = _change_univ_poly_ring(R, parent(p), cached)) where T
   upgrade!(p)
   return UniversalRingElem(change_base_ring(R, data(p); parent = base_ring(parent)), parent)
end

function change_coefficient_ring(R::Ring, p::UniversalRingElem{<:MPolyRingElem, T}; cached::Bool=true, parent::UniversalRing{<:MPolyRingElem} = _change_univ_poly_ring(R, parent(p), cached)) where T
  return change_base_ring(R, p, cached = cached, parent = parent)
end

################################################################################
#
#  Map
#
################################################################################

function map_coefficients(f::Any, p::UniversalRingElem{<:MPolyRingElem}; cached::Bool=true, parent::UniversalRing{<:MPolyRingElem} = _change_univ_poly_ring(parent(f(zero(coefficient_ring(p)))), parent(p), cached))
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
#   Unsafe functions
#
###############################################################################

function fit!(a::UniversalRingElem{<:MPolyRingElem}, n::Int)
  fit!(data(a), n)
end

###############################################################################
#
#   Parent object overload
#
###############################################################################

function (a::UniversalRing{<:MPolyRingElem, T})(b::Vector{T}, m::Vector{Vector{Int}}) where T
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

###############################################################################
#
#   constructors
#
###############################################################################

@doc raw"""
    universal_polynomial_ring(R::Ring, varnames::Vector{Symbol}; cached::Bool=true, internal_ordering::Symbol=:lex)
    universal_polynomial_ring(R::Ring; cached::Bool=true, internal_ordering::Symbol=:lex)

Given a coefficient ring `R` and variable names, say `varnames = [:x1, :x2, ...]`, return
a tuple `S, [x1, x2, ...]` of the universal polynomial ring `S = R[x1, x2, \dots]` and its generators `x1, x2, \dots`.

If `varnames` is omitted, return an object representing
the universal polynomial ring `S = R[\ldots]` with no variables in it initially.

# Examples

```jldoctest
julia> S, (x,y) = universal_polynomial_ring(ZZ, [:x,:y])
(Universal polynomial ring over Integers, UniversalRingElem{AbstractAlgebra.Generic.MPoly{BigInt}, BigInt}[x, y])

julia> z = gen(S, :z)
z

julia> x*y - z
x*y - z

julia> S = universal_polynomial_ring(ZZ)
Universal polynomial ring over Integers

julia> x = gen(S, :x)
x

julia> y, z = gens(S, [:y, :z])
2-element Vector{UniversalRingElem{AbstractAlgebra.Generic.MPoly{BigInt}, BigInt}}:
 y
 z

julia> x*y - z
x*y - z
```
"""
function universal_polynomial_ring(R::Ring, varnames::Vector{Symbol}; cached::Bool=true, internal_ordering::Symbol=:lex)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   S = get_cached!(UniversalPolyRingID, (R, internal_ordering), cached) do
      P = poly_ring(R, varnames; internal_ordering)
      T = elem_type(P)
      U = elem_type(coefficient_ring(P))
      UniversalRing{T, U}(P)
   end
   return (S, gens(S, varnames))
end

const UniversalPolyRingID = CacheDictType{Tuple{Ring, Symbol}, Ring}()

function universal_polynomial_ring(R::Ring; cached::Bool=true, internal_ordering::Symbol=:lex)
   return universal_polynomial_ring(R, Symbol[]; internal_ordering, cached)[1]
end

@varnames_interface universal_polynomial_ring(R::Ring, s)
