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

base_ring(S::UniversalPolyRing) = S.base_ring

coefficient_ring(S::UniversalPolyRing) = base_ring(S)

function is_domain_type(::Type{T}) where {S <: RingElement, U <: MPolyRingElem{S}, T <: UnivPoly{S, U}}
   return is_domain_type(S)
end

function is_exact_type(::Type{T}) where {S <: RingElement, U <: MPolyRingElem{S}, T <: UnivPoly{S, U}}
   return is_exact_type(S)
end

parent(p::UnivPoly) = p.parent

elem_type(::Type{UniversalPolyRing{T, U}}) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}} =
   UnivPoly{T, U}

parent_type(::Type{UnivPoly{T, U}}) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}} =
   UniversalPolyRing{T, U}

function mpoly_ring(S::UniversalPolyRing{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}}
   return S.mpoly_ring::Generic.MPolyRing{T}
end

number_of_variables(S::UniversalPolyRing) = length(S.S)

symbols(S::UniversalPolyRing) = S.S

function vars(p::UnivPoly{T, U}) where {T, U}
   S = parent(p)
   V = vars(p.p)
   return [UnivPoly{T, U}(v, S) for v in V]
end

ordering(p::UniversalPolyRing) = ordering(mpoly_ring(p))

function check_parent(a::UnivPoly{T, U}, b::UnivPoly{T, U}, throw::Bool = true) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}}
   flag = parent(a) != parent(b)
   flag & throw && error("Incompatible polynomial rings in polynomial operation")
   return !flag
end

###############################################################################
#
#   Manipulating terms and monomials
#
###############################################################################

exponent_vector(p::UnivPoly, i::Int) = exponent_vector(p.p, i)

function exponent(p::UnivPoly, i::Int, j::Int)
   if j > nvars(parent(p.p)) && j <= nvars(parent(p))
      return 0
   end
   return exponent(p.p, i, j)
end

function set_exponent_vector!(p::UnivPoly, i::Int, exps::Vector{Int})
   S = parent(p)
   len = length(exps)
   if len != nvars(parent(p.p))
      p.p = upgrade(S, p.p)
      if len < nvars(S)
         exps = vcat(exps, zeros(Int, nvars(S) - len))
      end
   end
   p.p = set_exponent_vector!(p.p, i, exps)
   return p
end

function coeff(p::UnivPoly, exps::Vector{Int})
   S = parent(p)
   len = length(exps)
   n = nvars(parent(p.p))
   if len > n
      if !iszero(exps[n + 1:len])
         return base_ring(S)()
      end
      return coeff(p.p, exps[1:n])
   end
   n = nvars(parent(p.p))
   if len < n
      exps = vcat(exps, zeros(Int, n - len))
   end
   return coeff(p.p, exps)
end

function setcoeff!(p::UnivPoly, exps::Vector{Int}, c::T) where T <: RingElement
   c = base_ring(p.p)(c)
   S = parent(p)
   len = length(exps)
   if len != nvars(parent(p.p))
      p.p = upgrade(S, p.p)
      if len < nvars(S)
         exps = vcat(exps, zeros(Int, nvars(S) - len))
      end
   end
   p.p = setcoeff!(p.p, exps, c)
   return p
end

function setcoeff!(a::UnivPoly{T, U}, i::Int, c::RingElement) where {T, U}
   setcoeff!(a.p, i, c)
   return a
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(R::UniversalPolyRing{T, U}) where {T, U} = UnivPoly{T, U}(zero(mpoly_ring(R)), R)

one(R::UniversalPolyRing{T, U}) where {T, U} = UnivPoly{T, U}(one(mpoly_ring(R)), R)

iszero(p::UnivPoly) = iszero(p.p)

isone(p::UnivPoly) = isone(p.p)

is_unit(p::UnivPoly) = is_unit(p.p)

is_zero_divisor(p::UnivPoly) = is_zero_divisor(p.p)

function is_zero_divisor_with_annihilator(p::UnivPoly{T, U}) where {T, U}
   f, b = is_zero_divisor_with_annihilator(p.p)
   return f, UnivPoly{T, U}(b, parent(p))
end

is_gen(p::UnivPoly) = is_gen(p.p)

is_homogeneous(p::UnivPoly) = is_homogeneous(p.p)

is_monomial(p::UnivPoly) = is_monomial(p.p)

is_constant(p::UnivPoly) = is_constant(p.p)

is_term(p::UnivPoly) = is_term(p.p)

coeff(p::UnivPoly, i::Int) = coeff(p.p, i)

function coeff(p::UnivPoly{T, U}, m::UnivPoly{T, U}) where {T, U}
   !is_monomial(m) && error("Not a monomial")
   check_parent(p, m)
   v1 = first(exponent_vectors(m))
   len = length(v1)
   n = nvars(parent(p.p))
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
   cvzip = zip(coefficients(p.p), exponent_vectors(p.p))
   for (c, v) in cvzip
      if v == v1
         return c
      end
   end
   return zero(R)
end

function monomial(p::UnivPoly{T, U}, i::Int) where {T, U}
   S = parent(p)
   m = monomial(p.p, i)
   return UnivPoly{T, U}(m, S)
end

function monomial!(m::UnivPoly{T, U}, p::UnivPoly{T, U}, i::Int) where {T, U}
   parent(m) != parent(p) && error("Incompatible monomial")
   if parent(m.p) != parent(p.p)
      m.p = parent(p.p)()
   end
   m.p = monomial!(m.p, p.p, i)
   return m
end

function term(p::UnivPoly{T, U}, i::Int) where {T, U}
   S = parent(p)
   t = term(p.p, i)
   return UnivPoly{T, U}(t, S)
end

leading_coefficient(p::UnivPoly) = leading_coefficient(p.p)

trailing_coefficient(p::UnivPoly) = trailing_coefficient(p.p)

function tail(p::UnivPoly{T, U}) where {T, U}
   S = parent(p)
   return UnivPoly{T, U}(tail(p.p), S)
end

constant_coefficient(p::UnivPoly) = constant_coefficient(p.p)

function leading_monomial(p::UnivPoly{T, U}) where {T, U}
   S = parent(p)
   return UnivPoly{T, U}(leading_monomial(p.p), S)
end

function leading_term(p::UnivPoly{T, U}) where {T, U}
   S = parent(p)
   return UnivPoly{T, U}(leading_term(p.p), S)
end

max_fields(p::UnivPoly) = max_fields(p.p)

function degree(p::UnivPoly, i::Int)
   if i <= nvars(parent(p)) && i > nvars(parent(p.p))
      return 0
   end
   return degree(p.p, i)
end

function degree(f::UnivPoly{T, U}, x::UnivPoly{T, U}) where {T, U}
   check_parent(f, x)
   return degree(f, var_index(x))
end

function degrees(p::UnivPoly)
   v = degrees(p.p)
   len = length(v)
   num = nvars(parent(p))
   if len < num
      v = vcat(v, zeros(Int, num - len))
   end
   return v
end

total_degree(p::UnivPoly) = total_degree(p.p)

length(p::UnivPoly) = length(p.p)

function gen(S::UniversalPolyRing{T, U}, s::VarName) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}}
   i = findfirst(==(Symbol(s)), S.S)
   if i === nothing
      push!(S.S, Symbol(s))
      S.mpoly_ring = MPolyRing{T}(S.base_ring, S.S, S.ord, false)
      i = length(S.S)
   end
   return UnivPoly{T, U}(gen(mpoly_ring(S), i), S)
end

gens(S::UniversalPolyRing, v::Vector{<:VarName}) = tuple([gen(S, s) for s in v]...)

function gen(S::UniversalPolyRing{T, U}, i::Int) where {T, U}
   @boundscheck 1 <= i <= nvars(S) || throw(ArgumentError("generator index out of range"))
   return UnivPoly{T, U}(gen(mpoly_ring(S), i), S)
end

function gens(S::UniversalPolyRing{T, U}) where {T, U}
   n = nvars(S)
   return UnivPoly{T, U}[gen(S, i) for i in 1:n]
end

var_index(x::UnivPoly) = var_index(x.p)

function vars(p::UnivPoly{T, U}) where {T <: RingElement, U}
   V = vars(p.p)
   S = parent(p)
   return [UnivPoly{T, U}(v, S) for v in V]
end

canonical_unit(p::UnivPoly) = canonical_unit(p.p)

characteristic(R::UniversalPolyRing) = characteristic(base_ring(R))

function Base.hash(p::UnivPoly, h::UInt)
   b = 0xcf418d4529109236%UInt
   for (c, v) in zip(coefficients(p.p), exponent_vectors(p.p))
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

function deepcopy_internal(p::UnivPoly{T, U}, dict::IdDict) where {T, U}
   return UnivPoly{T, U}(deepcopy_internal(p.p, dict), p.parent)
end

###############################################################################
#
#   Multivariate coefficients
#
###############################################################################

function coeff(p::UnivPoly{T, U}, vars::Vector{Int}, exps::Vector{Int}) where {T, U}
   len = length(vars)
   len != length(exps) && error("Number of variables does not match number of exponents")
   S = parent(p)
   n = nvars(S)
   num = nvars(parent(p.p))
   vars2 = Vector{Int}(undef, 0)
   exps2 = Vector{Int}(undef, 0)
   for i = 1:len
      vars[i] > n && error("Variable index not in range")
      if vars[i] <= num
         push!(vars2, vars[i])
         push!(exps2, exps[i])
      end
   end
   return UnivPoly{T, U}(coeff(p.p, vars2, exps2), S)
end

function coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where {U, V, T <: UnivPoly{U, V}}
   varidx = [var_index(x) for x in vars]
   return coeff(a, varidx, exps)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::UniversalPolyRing)
   print(io, "Universal Polynomial Ring over ")
   show(io, base_ring(R))
end

function Base.show(io::IO, ::MIME"text/plain", a::UnivPoly)
   print(io, AbstractAlgebra.obj_to_string(a.p, context = io))
end

function Base.show(io::IO, a::UnivPoly)
   print(io, AbstractAlgebra.obj_to_string(a.p, context = io))
end

###############################################################################
#
#   Unary operations
#
###############################################################################
   
function -(p::UnivPoly{T, U}) where {T, U}
   return UnivPoly{T, U}(-p.p, p.parent)
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function univ_promote(x::UnivPoly{T, U}, y::UnivPoly{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}}
   nx = nvars(parent(x.p))
   ny = nvars(parent(y.p))
   if nx == ny
      return x, y
   end
   S = parent(x)
   return S(x), S(y)
end

function +(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   a, b = univ_promote(a, b)
   return UnivPoly{T, U}(a.p + b.p, a.parent)
end

function -(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   a, b = univ_promote(a, b)
   return UnivPoly{T, U}(a.p - b.p, a.parent)
end

function *(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   a, b = univ_promote(a, b)
   return UnivPoly{T, U}(a.p*b.p, a.parent)
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

function Base.eltype(x::UnivPolyCoeffs{T}) where T <: AbstractAlgebra.UniversalPolyRingElem{S} where S <: RingElement
   return S
end

function Base.eltype(x::UnivPolyExponentVectors{T}) where T <: AbstractAlgebra.UniversalPolyRingElem{S} where S <: RingElement
   return Vector{Int}
end

function Base.eltype(x::UnivPolyMonomials{T}) where T <: AbstractAlgebra.UniversalPolyRingElem{S} where S <: RingElement
   return T
end

function Base.eltype(x::UnivPolyTerms{T}) where T <: AbstractAlgebra.UniversalPolyRingElem{S} where S <: RingElement
   return T
end

###############################################################################
#
#   Square root
#
###############################################################################

function Base.sqrt(p::UnivPoly{T, U}; check=true) where {T, U}
   S = parent(p)
   s = sqrt(p.p; check=check)
   return UnivPoly{T, U}(s, S)
end

function is_square(p::UnivPoly)
   return is_square(p.p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

function *(p::UnivPoly{T, U}, n::Union{Integer, Rational, AbstractFloat}) where {T, U}
   S = parent(p)
   return UnivPoly{T, U}(p.p*n, S)
end

function *(p::UnivPoly{T, U}, n::T) where {T <: RingElem, U}
   S = parent(p)
   return UnivPoly{T, U}(p.p*n, S)
end

*(n::Union{Integer, Rational, AbstractFloat}, p::UnivPoly) = p*n

*(n::T, p::UnivPoly{T, U}) where {T <: RingElem, U} = p*n

function divexact(p::UnivPoly{T, U}, n::Union{Integer, Rational, BigFloat}; check::Bool=true) where {T, U}
   S = parent(p)
   return UnivPoly{T, U}(divexact(p.p, n; check=check), S)
end

function divexact(p::UnivPoly{T, U}, n::T; check::Bool=true) where {T <: RingElem, U}
   S = parent(p)
   return UnivPoly{T, U}(divexact(p.p, n; check=check), S)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   if length(a) != length(b)
      return false
   end
   n1 = nvars(parent(a.p))
   n2 = nvars(parent(b.p))
   if n1 == n2
      for (v1, v2) in zip(exponent_vectors(a.p), exponent_vectors(b.p))
         if v1 != v2
            return false
         end
      end
   elseif n1 > n2
      for (v1, v2) in zip(exponent_vectors(a.p), exponent_vectors(b.p))
         if v1[1:n2] != v2 || !iszero(v1[n2 + 1:n1])
            return false
         end
      end
   else # n2 > n1
      for (v1, v2) in zip(exponent_vectors(a.p), exponent_vectors(b.p))
         if v1 != v2[1:n1] || !iszero(v2[n1 + 1:n2])
            return false
         end
      end
   end
   for (c1, c2) in zip(coefficients(a.p), coefficients(b.p))
      if c1 != c2
         return false
      end
   end
   return true
end

function isless(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   if parent(a.p) === parent(b.p)
      return isless(a.p, b.p)
   end
   S = parent(a)
   num = nvars(S)
   s = a.p
   t = b.p
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

==(p::UnivPoly, n::Union{Integer, Rational, AbstractFloat}) = p.p == n

==(n::Union{Integer, Rational, AbstractFloat}, p::UnivPoly) = p.p == n

==(p::UnivPoly{T, U}, n::T) where {T <: RingElem, U} = p.p == n

==(n::T, p::UnivPoly{T, U}) where {T <: RingElem, U} = p.p == n

###############################################################################
#
#   Powering
#
###############################################################################

function ^(p::UnivPoly{T, U}, b::Int) where {T, U}
   S = parent(p)
   return UnivPoly{T, U}(p.p^b, S)
end

###############################################################################
#
#   Inflation/deflation
#
###############################################################################

deflation(p::UnivPoly{T, U}) where {T, U} = deflation(p.p)

function deflate(p::UnivPoly{T, U}, shift::Vector{Int}, defl::Vector{Int}) where {T, U}
   S = parent(p)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   num = nvars(parent(p.p))
   pp = p.p
   if vlen == num
      return UnivPoly{T, U}(deflate(pp, shift, defl), S)
   end
   if vlen > num
      pp = upgrade(S, pp)
      num = nvars(parent(pp))
   end
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return UnivPoly{T, U}(deflate(pp, shift, defl), S)
end

function deflate(p::UnivPoly{T, U}, defl::Vector{Int}) where {T, U}
   return deflate(p, zeros(Int, length(defl)), defl)
end

function inflate(p::UnivPoly{T, U}, shift::Vector{Int}, defl::Vector{Int}) where {T, U}
   S = parent(p)
   vlen = length(shift)
   vlen != length(defl) && error("Vector lengths do not match")
   num = nvars(parent(p.p))
   pp = p.p
   if vlen == num
      return UnivPoly{T, U}(inflate(pp, shift, defl), S)
   end
   if vlen > num
      pp = upgrade(S, pp)
      num = nvars(parent(pp))
   end
   if vlen < num
      shift = vcat(shift, zeros(Int, num - vlen))
      defl = vcat(defl, ones(Int, num - vlen))
   end
   return UnivPoly{T, U}(inflate(pp, shift, defl), S)
end

function inflate(p::UnivPoly{T, U}, defl::Vector{Int}) where {T, U}
   return inflate(p, zeros(Int, length(defl)), defl)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::UnivPoly{T, U}, b::UnivPoly{T, U}; check::Bool=true) where {T, U}
   check_parent(a, b)
   a, b = univ_promote(a, b)
   return UnivPoly{T, U}(divexact(a.p, b.p; check=check), a.parent)
end

function divides(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   a, b = univ_promote(a, b)
   flag, q = divides(a.p, b.p)
   return flag, UnivPoly{T, U}(q, a.parent)
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function Base.div(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   a, b = univ_promote(a, b)
   return UnivPoly{T, U}(div(a.p, b.p), a.parent)
end

function Base.divrem(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   a, b = univ_promote(a, b)
   q, r = divrem(a.p, b.p)
   return UnivPoly{T, U}(q, a.parent), UnivPoly{T, U}(r, a.parent)
end

###############################################################################
#
#   Derivative
#
###############################################################################

function derivative(p::UnivPoly{T, U}, j::Int) where {T, U}
   j > nvars(parent(p)) && error("No such variable")
   if j > nvars(parent(p.p))
      return zero(parent(p))
   end
   return UnivPoly{T, U}(derivative(p.p, j), p.parent)
end

function derivative(p::UnivPoly{T, U}, x::UnivPoly{T, U}) where {T, U}
   return derivative(p, var_index(x))
end

###############################################################################
#
#   Remove and valuation
#
###############################################################################

function remove(z::UnivPoly{T, U}, p::UnivPoly{T, U}) where {T, U}
   check_parent(z, p)
   S = parent(z)
   z, p = univ_promote(z, p)
   val, q = remove(z.p, p.p)
   return val, UnivPoly{T, U}(q, S)
end

function valuation(z::UnivPoly{T, U}, p::UnivPoly{T, U}) where {T, U}
  v, _ = remove(z, p)
  return v
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(a::UnivPoly{T, U}, A::Vector{T}) where {T <: RingElem, U}
   R = base_ring(a)
   n = length(A)
   num = nvars(parent(a.p))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
      if nvars(parent(a.p)) == 0
         return constant_coefficient(a.p)*one(parent(A[1]))
      end
      return evaluate(a.p, A[1:num])
   end
   if n < num
      A = vcat(A, [zero(R) for i = 1:num - n])
   end
   return evaluate(a.p, A)
end

function evaluate(a::UnivPoly{T, U}, A::Vector{V}) where {T <: RingElement, U, V <: Union{Integer, Rational, AbstractFloat}}
   n = length(A)
   num = nvars(parent(a.p))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
      if nvars(parent(a.p)) == 0
         return constant_coefficient(a.p)*one(parent(A[1]))
      end
      return evaluate(a.p, A[1:num])
   end
   if n < num
      A = vcat(A, zeros(V, num - n))
   end
   return evaluate(a.p, A)
end

function evaluate(a::UnivPoly{T, U}, A::Vector{V}) where {T <: RingElement, U, V <: RingElement}
   n = length(A)
   num = nvars(parent(a.p))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
      if nvars(parent(a.p)) == 0
         return constant_coefficient(a.p)*one(parent(A[1]))
      end
      return evaluate(a.p, A[1:num])
   end
   if n < num
      if n == 0
         R = base_ring(a)
         return evaluate(a.p, zeros(R, num))
      else
         R = parent(A[1])
         A = vcat(A, zeros(R, num - n))
         return evaluate(a.p, A)
     end
   end
   return evaluate(a.p, A)
end

function (a::UnivPoly{T, U})() where {T <: RingElement, U}
   return evaluate(a, T[])
end

function (a::UnivPoly{T, U})(vals::T...) where {T <: RingElement, U}
   return evaluate(a, [vals...])
end

function (a::UnivPoly{T, U})(val::V, vals::V...) where {T <: RingElement, U, V <: Union{Integer, Rational, AbstractFloat}}
   return evaluate(a, [val, vals...])
end

function (a::UnivPoly{T, U})(vals::Union{NCRingElem, RingElement}...) where {T <: RingElement, U}
   A = [vals...]
   n = length(vals)
   num = nvars(parent(a.p))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
      if nvars(parent(a.p)) == 0
         return constant_coefficient(a.p)*one(parent(A[1]))
      end
      return a.p(vals[1:num]...)
   end
   if n < num
      A = vcat(A, zeros(Int, num - n))
   end
   return a.p(A...)
end

function evaluate(a::UnivPoly{T, U}, vals::Vector{V}) where {T <: RingElement, U, V <: NCRingElem}
   return a(vals...)
end

function evaluate(a::UnivPoly{T, U}, vars::Vector{Int}, vals::Vector{V}) where {T <: RingElement, U, V <: RingElement}
   length(vars) != length(vals) && error("Numbers of variables and values do not match")
   vars2 = Vector{Int}(undef, 0)
   vals2 = Vector{U}(undef, 0)
   num = nvars(parent(a.p))
   S = parent(a)
   n = nvars(S)
   for i = 1:length(vars)
      vars[i] > n && error("Unknown variable")
      if vars[i] <= num
         push!(vars2, vars[i])
         push!(vals2, S(vals[i]).p)
      end
   end
   return UnivPoly(evaluate(S(a).p, vars2, vals2), S)
end

function evaluate(a::S, vars::Vector{S}, vals::Vector{V}) where {S <: UnivPoly{T, U}, V <: RingElement} where {T <: RingElement, U}
   varidx = Int[var_index(x) for x in vars]
   return evaluate(a, varidx, vals)
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T <: RingElement, U}
   check_parent(a, b)
   a, b = univ_promote(a, b)
   return UnivPoly{T, U}(gcd(a.p, b.p), a.parent)
end

function lcm(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T <: RingElement, U}
   check_parent(a, b)
   a, b = univ_promote(a, b)
   return UnivPoly{T, U}(lcm(a.p, b.p), a.parent)
end

###############################################################################
#
#   Univariate polynomials
#
###############################################################################

function to_univariate(R::AbstractAlgebra.PolyRing{T}, p::UnivPoly{T, U}) where {T <: RingElement, U}
   return to_univariate(R, p.p)
end

is_univariate(p::UnivPoly) = is_univariate(p.p)

is_univariate(R::UniversalPolyRing) = is_univariate(mpoly_ring(R))

function coefficients_of_univariate(p::UnivPoly, check_univariate::Bool=true)
   return coefficients_of_univariate(p.p, check_univariate)
end

################################################################################
#
#  Change base ring
#
################################################################################

function _change_univ_poly_ring(R, Rx, cached::Bool)
   P, _ = AbstractAlgebra.polynomial_ring(R, map(string, symbols(Rx)), ordering = ordering(Rx), cached = cached)
   S = AbstractAlgebra.UniversalPolynomialRing(R; ordering=ordering(Rx), cached=cached)
   S.S = deepcopy(symbols(Rx))
   S.mpoly_ring = P
   return S
end

function change_base_ring(R::Ring, p::UnivPoly{T, U}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(R, parent(p), cached)) where {T <: RingElement, U}
   base_ring(parent) != R && error("Base rings do not match.")
   return _map(R, p, parent)
end

function change_coefficient_ring(R::Ring, p::UnivPoly{T, U}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(R, parent(p), cached)) where {T <: RingElement, U}
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

function sort_terms!(p::UnivPoly{T, U}) where {T, U}
   p.p = sort_terms!(p.p)
   return p
end

function combine_like_terms!(p::UnivPoly{T, U}) where {T, U}
   p.p = combine_like_terms!(p.p)
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
   rand(GLOBAL_RNG, S, term_range, exp_bound, v...)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!(a::UnivPoly, n::Int)
   fit!(a.p, n)
end

function add!(a::UnivPoly{T, U}, b::UnivPoly{T, U}, c::UnivPoly{T, U}) where {T <: RingElement, U}
   a.p = (b + c).p
   return a
end

function mul!(a::UnivPoly{T, U}, b::UnivPoly{T, U}, c::UnivPoly{T, U}) where {T <: RingElement, U}
   a.p = (b*c).p
   return a
end

function addeq!(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T <: RingElement, U}
   a.p = (a + b).p
   return a
end

function addmul!(a::UnivPoly{T, U}, b::UnivPoly{T, U}, c::UnivPoly{T, U}) where {T <: RingElement, U}
   a.p = (a + b*c).p
   return a
end

function zero!(a::UnivPoly{T, U}) where {T <: RingElement, U}
   a.p = zero!(a.p)
   return a
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{UnivPoly{T, U}}, ::Type{UnivPoly{T, U}}) where {T <: RingElement, U} = UnivPoly{T, U}

function promote_rule(::Type{UnivPoly{T, U}}, ::Type{V}) where {T <: RingElement, U, V <: RingElement}
   promote_rule(T, V) == T ? UnivPoly{T, U} : Union{}
end

###############################################################################
#
#   Parent object overload
#
###############################################################################

function upgrade(S::UniversalPolyRing{T, U}, pp::MPoly{T}) where {T, U}
   alloc = length(pp.coeffs)
   n = nvars(S) - nvars(parent(pp))
   ctx = MPolyBuildCtx(mpoly_ring(S))
   v0 = zeros(Int, n)
   for (c, v) in zip(coefficients(pp), exponent_vectors(pp))
      push_term!(ctx, c, vcat(v, v0))
   end
   p = finish(ctx)
   fit!(p, alloc)
   return p
end

function upgrade(S::UniversalPolyRing{T, U}, pp::MPolyRingElem{T}) where {T, U}
   n = nvars(S) - nvars(parent(pp))
   ctx = MPolyBuildCtx(mpoly_ring(S))
   v0 = zeros(Int, n)
   for (c, v) in zip(coefficients(pp), exponent_vectors(pp))
      push_term!(ctx, c, vcat(v, v0))
   end
   return finish(ctx)
end

function (a::UniversalPolyRing{T, U})(b::RingElement) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}}
   return a(base_ring(a)(b))
end

function (a::UniversalPolyRing{T, U})() where {T <: RingElement, U}
   return UnivPoly{T, U}(mpoly_ring(a)(), a)
end

function (a::UniversalPolyRing{T, U})(b::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}}
   return UnivPoly{T, U}(mpoly_ring(a)(b), a)
end

function (a::UniversalPolyRing{T, U})(b::T) where {T <: RingElem, U <: AbstractAlgebra.MPolyRingElem{T}}
   return UnivPoly{T, U}(mpoly_ring(a)(b), a)
end

function (S::UniversalPolyRing{T, U})(p::UnivPoly{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}}
   parent(p) !== S && error("Unable to coerce")
   n = nvars(S) - nvars(parent(p.p))
   if n != 0
      p = UnivPoly{T, U}(upgrade(S, p.p), S)
   end
   return p
end

function (a::UniversalPolyRing{T, U})(b::Vector{T}, m::Vector{Vector{Int}}) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}}
   if length(m) != 0
      len = length(m[1])
      num = nvars(mpoly_ring(a))
      if len != num
         for i = 1:length(m)
            m[i] = vcat(m[i], zeros(Int, num - len))
         end
      end
   end
   return UnivPoly{T, U}(mpoly_ring(a)(b, m), a)
end

###############################################################################
#
#   UniversalPolynomialRing constructor
#
###############################################################################

function UniversalPolynomialRing(R::Ring; ordering=:lex, cached::Bool=true)
   T = elem_type(R)
   U = Generic.MPoly{T}

   return UniversalPolyRing{T, U}(R, ordering, cached)
end

