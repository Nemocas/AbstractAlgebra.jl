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

base_ring(S::UnivPolyRing) = S.base_ring

parent(p::UnivPoly) = p.parent

elem_type(::Type{UnivPolyRing{T, U}}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}} =
   UnivPoly{T, U}

parent_type(::Type{UnivPoly{T, U}}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}} =
   UnivPolyRing{T, U}

function mpoly_ring(S::UnivPolyRing{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
   return S.mpoly_ring::parent_type(mpoly_type(T))
end

nvars(S::UnivPolyRing) = length(S.S)

symbols(S::UnivPolyRing) = S.S

function vars(p::UnivPoly{T, U}) where {T, U}
   S = parent(p)
   V = vars(p.p)
   return [UnivPoly{T, U}(v, S) for v in V]
end

ordering(p::UnivPolyRing) = ordering(mpoly_ring(p))

function check_parent(a::UnivPoly{T, U}, b::UnivPoly{T, U}, throw::Bool = true) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
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
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

isunit(p::UnivPoly) = isunit(p.p)

isgen(p::UnivPoly) = isgen(p.p)

ishomogeneous(p::UnivPoly) = ishomogeneous(p.p)

coeff(p::UnivPoly, i::Int) = coeff(p.p, i)

trailing_coefficient(p) = trailing_coefficient(p.p)

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

max_fields(p::UnivPoly) = max_fields(p.p)

function degree(p::UnivPoly, i::Int)
   if i <= nvars(parent(p)) && i > nvars(parent(p.p))
      return 0
   end
   return degree(p.p, i)
end

total_degree(p::UnivPoly) = total_degree(p.p)

length(p::UnivPoly) = length(p.p)

function gen(S::UnivPolyRing{T, U}, s::Symbol) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
   i = findfirst(x->x==s, S.S)
   if typeof(i) == Nothing
      push!(S.S, s)
      S.mpoly_ring = PolynomialRing(base_ring(S), S.S; cached=true, ordering=S.ord)[1]
      i = length(S.S)
   end
   return UnivPoly{T, U}(gen(mpoly_ring(S), i), S)
end

gen(S::UnivPolyRing, s::Union{Char, String}) = gen(S, Symbol(s))

gens(S::UnivPolyRing, v::Vector{Symbol}) = tuple([gen(S, s) for s in v]...)

gens(S::UnivPolyRing, v::Vector{T}) where T <: Union{Char, String} = gens(S, [Symbol(s) for s in v])

var_index(x::UnivPoly) = var_index(x.p)

function Base.hash(p::UnivPoly, h::UInt)
   b = 0xcf418d4529109236%UInt
   for (c, v) in zip(coefficients(p.p), exponent_vectors(p.p))
      b = xor(b, xor(Base.hash(v[1:findlast(x->x>0, v)], h), h))
      b = xor(b, xor(hash(c, h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, p::UnivPolyRing)
   print(io, "Universal Polynomial Ring")
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

function univ_promote(x::UnivPoly{T, U}, y::UnivPoly{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
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

function Base.eltype(x::UnivPolyCoeffs{T}) where T <: AbstractAlgebra.UnivPolyElem{S} where S <: RingElement
   return S
end

function Base.eltype(x::UnivPolyExponentVectors{T}) where T <: AbstractAlgebra.UnivPolyElem{S} where S <: RingElement
   return Vector{Int}
end

function Base.eltype(x::UnivPolyMonomials{T}) where T <: AbstractAlgebra.UnivPolyElem{S} where S <: RingElement
   return T
end

function Base.eltype(x::UnivPolyTerms{T}) where T <: AbstractAlgebra.UnivPolyElem{S} where S <: RingElement
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

function issquare(p::UnivPoly)
   return issquare(p.p)
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
#   Evaluation
#
###############################################################################

function evaluate(a::UnivPoly{T, U}, A::Vector{T}) where {T <: RingElement, U}
   R = base_ring(a)
   n = length(A)
   num = nvars(parent(a.p))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
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

function (a::UnivPoly{T, U})(vals::T...) where {T <: RingElement, U}
   return evaluate(a, [vals...])
end

function (a::UnivPoly{T, U})(vals::V...) where {T <: RingElement, U, V <: Union{Integer, Rational, AbstractFloat}}
   return evaluate(a, [vals...])
end

function (a::UnivPoly{T, U})(vals::Union{NCRingElem, RingElement}...) where {T <: RingElement, U}
   A = [vals...]
   n = length(vals)
   num = nvars(parent(a.p))
   if n > num
      n > nvars(parent(a)) && error("Too many values")
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

function evaluate(a::UnivPoly{T}, vars::Vector{Int}, vals::Vector{U}) where {T <: RingElement, U <: RingElement}
   length(vars) != length(vals) && error("Numbers of variables and values do not match")
   vars2 = Vector{Int}(undef, 0)
   vals2 = Vector{U}(undef, 0)
   num = nvars(parent(a.p))
   n = nvars(parent(a))
   for i = 1:length(vars)
      vars[i] > n && error("Unknown variable")
      if vars[i] <= num
         push!(vars2, vars[i])
         push!(vals2, vals[i])
      end
   end
   return evaluate(a.p, vars2, vals2)
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
#   Unsafe functions
#
###############################################################################

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

function upgrade(S::UnivPolyRing{T, U}, pp::MPolyElem{T}) where {T, U}
   n = nvars(S) - nvars(parent(pp))
   ctx = MPolyBuildCtx(mpoly_ring(S))
   v0 = zeros(Int, n)
   for (c, v) in zip(coefficients(pp), exponent_vectors(pp))
      push_term!(ctx, c, vcat(v, v0))
   end
   return finish(ctx)
end

function (a::UnivPolyRing{T, U})() where {T <: RingElement, U}
   return UnivPoly{T, U}(mpoly_ring(a)(), a)
end

function (S::UnivPolyRing{T, U})(p::UnivPoly{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
   parent(p) !== S && error("Unable to coerce")
   n = nvars(S) - nvars(parent(p.p))
   if n != 0
      p = UnivPoly{T, U}(upgrade(S, p.p), S)
   end
   return p
end

###############################################################################
#
#   UniversalPolynomialRing constructor
#
###############################################################################

function UniversalPolynomialRing(R::Ring; ordering=:lex, cached=true)
   T = elem_type(R)
   U = mpoly_type(T)

   return UnivPolyRing{T, U}(R, ordering)
end

