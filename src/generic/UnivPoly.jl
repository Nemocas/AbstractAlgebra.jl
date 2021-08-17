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

