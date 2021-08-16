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

###############################################################################
#
#   Basic manipulation
#
###############################################################################

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

function check_parent(a::UnivPoly{T, U}, b::UnivPoly{T, U}, throw::Bool = true) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
   flag = parent(a) != parent(b)
   flag & throw && error("Incompatible polynomial rings in polynomial operation")
   return !flag
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
#   Parent object overload
#
###############################################################################

function (S::UnivPolyRing{T, U})(p::UnivPoly{T, U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
   parent(p) !== S && error("Unable to coerce")
   n = nvars(S) - nvars(parent(p.p))
   if n != 0
      ctx = MPolyBuildCtx(mpoly_ring(S))
      v0 = zeros(Int, n)
      for (c, v) in zip(coefficients(p.p), exponent_vectors(p.p))
         push_term!(ctx, c, vcat(v, v0))
      end
      p = UnivPoly{T, U}(finish(ctx), S)
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

