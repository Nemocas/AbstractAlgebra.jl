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

base_ring_type(::Type{<:UniversalPolyRing{T}}) where T = mpoly_ring_type(T)
base_ring(S::UniversalPolyRing) = S.base_ring::base_ring_type(S)

coefficient_ring_type(::Type{<:UniversalPolyRing{T}}) where T = parent_type(T)
coefficient_ring(S::UniversalPolyRing) = coefficient_ring(base_ring(S))::coefficient_ring_type(S)

function is_domain_type(::Type{<:UnivPoly{S}}) where {S <: RingElement}
   return is_domain_type(S)
end

function is_exact_type(::Type{<:UnivPoly{S}}) where {S <: RingElement}
   return is_exact_type(S)
end

parent(p::UnivPoly) = p.parent::parent_type(p)

elem_type(::Type{UniversalPolyRing{T}}) where {T<:RingElement} = UnivPoly{T}

parent_type(::Type{UnivPoly{T}}) where {T<:RingElement} = UniversalPolyRing{T}

number_of_variables(S::UniversalPolyRing) = number_of_variables(base_ring(S))

number_of_generators(S::UniversalPolyRing) = number_of_generators(base_ring(S))

symbols(S::UniversalPolyRing) = symbols(base_ring(S))

internal_ordering(p::UniversalPolyRing) = internal_ordering(base_ring(p))

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
   upgrade!(p)
   if len < nvars(S)
      exps = vcat(exps, zeros(Int, nvars(S) - len))
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

function setcoeff!(p::UnivPoly, exps::Vector{Int}, c::T) where T <: RingElement
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

function setcoeff!(a::UnivPoly{T}, i::Int, c::RingElement) where {T}
   setcoeff!(data(a), i, c)
   return a
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(R::UniversalPolyRing{T}) where {T} = UnivPoly{T}(zero(base_ring(R)), R)

one(R::UniversalPolyRing{T}) where {T} = UnivPoly{T}(one(base_ring(R)), R)

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
   p, m = univ_promote(p, m)
   return coeff(data(p), data(m))
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
      S.base_ring = AbstractAlgebra._add_gens(base_ring(S), added_symbols)
   end
   return idx
end

function gen(S::UniversalPolyRing, s::VarName)
   i = findfirst(==(Symbol(s)), symbols(S))
   if i === nothing
      S.base_ring = AbstractAlgebra._add_gens(base_ring(S), [Symbol(s)])
      i = length(symbols(S))
   end
   return @inbounds gen(S, i)
end

function gen(S::UniversalPolyRing{T}, i::Int) where {T}
   @boundscheck 1 <= i <= nvars(S) || throw(ArgumentError("generator index out of range"))
   return UnivPoly{T}(gen(base_ring(S), i), S)
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
   return nothing, [UnivPoly{T}(gen(base_ring(S), i), S) for i in idx]
end

AbstractAlgebra.@varnames_interface _univ_poly_gens(R::UniversalPolyRing{T}, s) where {T}

function gens(S::UniversalPolyRing{T}, a::AbstractAlgebra.VarNames, as::AbstractAlgebra.VarNames...) where {T}
   res = _univ_poly_gens(S, a, as...)
   length(res) == 2 && return res[2] # special case for improved backwards compatibility
   return res[2:end]
end

var_index(x::UnivPoly) = var_index(data(x))

var_indices(p::UnivPoly) = var_indices(data(p))

function vars(p::UnivPoly{T}) where {T <: RingElement}
   V = vars(data(p))
   S = parent(p)
   return [UnivPoly{T}(v, S) for v in V]
end

canonical_unit(p::UnivPoly) = canonical_unit(data(p))

characteristic(R::UniversalPolyRing) = characteristic(coefficient_ring(R))

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
   show(io, coefficient_ring(R))
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
    function $op(p::UnivPoly{T}, n::JuliaRingElement) where {T}
       S = parent(p)
       return UnivPoly{T}($op(data(p),n), S)
    end

    function $op(p::UnivPoly{T}, n::T) where {T <: RingElem}
       S = parent(p)
       return UnivPoly{T}($op(data(p),n), S)
    end

    function $op(n::JuliaRingElement, p::UnivPoly{T}) where {T}
       S = parent(p)
       return UnivPoly{T}($op(n,data(p)), S)
    end

    function $op(n::T, p::UnivPoly{T}) where {T <: RingElem}
       S = parent(p)
       return UnivPoly{T}($op(n,data(p)), S)
    end
  end
end

function divexact(p::UnivPoly{T}, n::JuliaRingElement; check::Bool=true) where {T}
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
   return isless(data(upgrade!(a)), data(upgrade!(b)))
end

###############################################################################
#
#   Ad hoc comparison functions
#
###############################################################################

==(p::UnivPoly, n::JuliaRingElement) = data(p) == n

==(n::JuliaRingElement, p::UnivPoly) = data(p) == n

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
   num = nvars(S)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   upgrade!(p)
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return UnivPoly{T}(deflate(data(p), shift, defl), S)
end

function deflate(p::UnivPoly{T}, defl::Vector{Int}) where {T}
   return deflate(p, zeros(Int, length(defl)), defl)
end

function inflate(p::UnivPoly{T}, shift::Vector{Int}, defl::Vector{Int}) where {T}
   S = parent(p)
   num = nvars(S)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   upgrade!(p)
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return UnivPoly{T}(inflate(data(p), shift, defl), S)
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

function evaluate(a::UnivPoly, A::Vector{<:NCRingElement})
   isempty(A) && error("Evaluating at an empty list of values is not allowed")
   a2 = data(a)
   varidx = var_indices(a2)
   isempty(varidx) && return constant_coefficient(a2) # TODO: this is weird
   vals = zeros(parent(A[1]), nvars(parent(a2)))
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

function evaluate(a::UnivPoly{T}, vars::Vector{Int}, vals::Vector{V}) where {T <: RingElement, V <: RingElement}
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
   return UnivPoly(evaluate(a2, vars2, vals2), S)
end

function evaluate(a::S, vars::Vector{S}, vals::Vector{V}) where {S <: UnivPoly{T}, V <: RingElement} where {T <: RingElement}
   varidx = Int[var_index(x) for x in vars]
   return evaluate(a, varidx, vals)
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

is_univariate(R::UniversalPolyRing) = is_univariate(base_ring(R))

function coefficients_of_univariate(p::UnivPoly, check_univariate::Bool=true)
   return coefficients_of_univariate(data(p), check_univariate)
end

################################################################################
#
#  Change base ring
#
################################################################################

_change_univ_poly_ring(R, Rx, cached::Bool) = universal_polynomial_ring(R, symbols(Rx); internal_ordering=internal_ordering(Rx), cached)[1]

function change_base_ring(R::Ring, p::UnivPoly{T}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(R, parent(p), cached)) where {T <: RingElement}
   upgrade!(p)
   return UnivPoly(change_base_ring(R, data(p); parent = base_ring(parent)), parent)
end

function change_coefficient_ring(R::Ring, p::UnivPoly{T}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(R, parent(p), cached)) where {T <: RingElement}
  return change_base_ring(R, p, cached = cached, parent = parent)
end

################################################################################
#
#  Map
#
################################################################################

function map_coefficients(f::T, p::UnivPoly; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(parent(f(zero(coefficient_ring(p)))), parent(p), cached)) where T
   upgrade!(p)
   return UnivPoly(map_coefficients(f, data(p); parent = base_ring(parent)), parent)
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
   R = coefficient_ring(S)
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
   R = coefficient_ring(S)
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

function upgrade!(p::UnivPoly)
   R = base_ring(parent(p))
   if R != parent(p.p)
      p.p = AbstractAlgebra._upgrade(p.p, R)
   end
   return p
end

function (a::UniversalPolyRing{T})(b::RingElement) where {T <: RingElement}
   return a(coefficient_ring(a)(b))
end

function (a::UniversalPolyRing{T})() where {T <: RingElement}
   return UnivPoly{T}(base_ring(a)(), a)
end

function (a::UniversalPolyRing{T})(b::JuliaRingElement) where {T <: RingElement}
   return UnivPoly{T}(base_ring(a)(b), a)
end

function (a::UniversalPolyRing{T})(b::T) where {T <: RingElem}
   return UnivPoly{T}(base_ring(a)(b), a)
end

function (S::UniversalPolyRing{T})(p::UnivPoly{T}) where {T <: RingElement}
   parent(p) !== S && error("Unable to coerce")
   return upgrade!(p)
end

function (a::UniversalPolyRing{T})(b::Vector{T}, m::Vector{Vector{Int}}) where {T <: RingElement}
   if length(m) != 0
      len = length(m[1])
      num = nvars(base_ring(a))
      if len != num
         for i = 1:length(m)
            m[i] = vcat(m[i], zeros(Int, num - len))
         end
      end
   end
   return UnivPoly{T}(base_ring(a)(b, m), a)
end
