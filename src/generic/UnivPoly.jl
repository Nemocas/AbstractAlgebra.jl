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

elem_type(S::UnivPolyRing{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}} =
   UnivPoly{T, U}

parent_type(p::UnivPoly{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}} =
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

function promote(x::UnivPoly{T, U}, y::UnivPoly{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
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
   a, b = promote(a, b)
   return UnivPoly{T, U}(a.p + b.p, a.parent)
end

function -(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   a, b = promote(a, b)
   return UnivPoly{T, U}(a.p - b.p, a.parent)
end

function *(a::UnivPoly{T, U}, b::UnivPoly{T, U}) where {T, U}
   check_parent(a, b)
   a, b = promote(a, b)
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

