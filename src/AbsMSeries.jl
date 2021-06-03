###############################################################################
#
#   AbsMSeries.jl : Multivariate power series over rings, capped absolute
#                   precision
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

function O(a::AbsMSeriesElem{T}) where T <: RingElement
    if iszero(a)
       return deepcopy(a)
    end
    R = parent(a)
    p = poly(a)
    v = vars(p)
    (length(v) != 1 || length(p) != 1 || !isone(leading_coefficient(p))) &&
                                               error("Not a pure power in O()")
    ind = var_index(v[1])
    exps = first(exponent_vectors(p))
    prec = [i == ind ? exps[i] : R.prec_max[i] for i in 1:length(exps)]
    return R(parent(p)(), prec)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc Markdown.doc"""
    symbols(R::MSeriesRing)

Return a vector of symbols, one for each of the variables of the series ring
$R$.
"""
symbols(R::MSeriesRing) = R.sym

parent(a::MSeriesElem) = a.parent

function base_ring(R::MSeriesRing{T}) where T <: RingElement
    return base_ring(poly_ring(R))::parent_type(T)
end

base_ring(a::MSeriesElem) = base_ring(parent(a))

@doc Markdown.doc"""
    characteristic(a::MSeriesRing)

Return the characteristic of the base ring of the series `a`.
"""
function characteristic(a::MSeriesRing)
    return characteristic(base_ring(a))
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(a::AbsMSeriesElem,
                                     x = symbols(parent(a)); context = nothing)
   apoly = poly(a)

   poly_sum = Expr(:call, :+)
   n = nvars(parent(apoly))

   iter = zip(coefficients(apoly), exponent_vectors(apoly))
   cv = reverse!(collect(iter))

   for (c, v) in cv
      prod = Expr(:call, :*)
      if !isone(c)
         push!(prod.args, expressify(c, context = context))
      end
      for i in n:-1:1
         if v[i] > 1
            push!(prod.args, Expr(:call, :^, x[i], v[i]))
         elseif v[i] == 1
            push!(prod.args, x[i])
         end
      end
      push!(poly_sum.args, prod)
   end

   sum = Expr(:call, :+)

   push!(sum.args, poly_sum)

   for i in nvars(parent(a)):-1:1
      push!(sum.args, Expr(:call, :O, Expr(:call, :^, x[i], a.prec[i])))
   end

   return sum
end

function Base.show(io::IO, a::MSeriesElem)
   print(io, obj_to_string(a, context = io))
end

function Base.show(io::IO, ::MIME"text/plain", a::MSeriesElem)
   print(io, obj_to_string(a, context = io))
end

function show(io::IO, a::MSeriesRing)
   v = join([String(s) for s in symbols(a)], ", ")
   print(io, "Multivariate power series ring in ", v, " over ")
   print(IOContext(io, :compact => true), base_ring(a))
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::MSeriesRing, _, _) = elem_type(S)

function RandomExtensions.make(S::MSeriesRing,
                                      term_range::UnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, vs[1])
   else
      make(S, term_range, make(R, vs...))
   end
end

function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make3{
                  <:RingElement,<:MSeriesRing,UnitRange{Int}}})
   S, term_range, v = sp[][1:end]
   f = S()
   g = gens(S)
   R = base_ring(S)
   prec = max_precision(S)
   for i = 1:rand(rng, term_range)
      term = S(1)
      for j = 1:length(g)
         term *= g[j]^rand(rng, 0:prec[j])
      end
      term *= rand(rng, v)
      f += term
   end
   return f
end

function rand(rng::AbstractRNG, S::MSeriesRing,
                                             term_range::UnitRange{Int}, v...)
   rand(rng, make(S, term_range, v...))
end

@doc Markdown.doc"""
    rand(S::MSeriesRing, term_range, v...)

Return a random element of the series ring $S$ with number of terms in the
range given by `term_range` and where coefficients of the series are randomly
generated in the base ring using the data given by `v`. The exponents of the
variable in the terms will be less than the precision caps for the Ring $S$
when it was created.
"""
function rand(S::MSeriesRing, term_range, v...)
   rand(GLOBAL_RNG, S, term_range, v...)
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

function PowerSeriesRing(R::Ring, prec::Vector{Int},
                  s::Vector{T}; cached=true, model=:capped_absolute) where
                                                                    T <: Symbol
   return Generic.PowerSeriesRing(R, prec, s; cached=cached, model=model)
end

function PowerSeriesRing(R::Ring, prec::Vector{Int},
   s::Vector{T}; cached=true, model=:capped_absolute) where
                                               T <: Union{Char, AbstractString}
   sym = [Symbol(v) for v in s]
   return Generic.PowerSeriesRing(R, prec, sym; cached=cached, model=model)
end

function PowerSeriesRing(R::Ring, prec::Int,
                  s::Vector{T}; cached=true, model=:capped_absolute) where
                                                                    T <: Symbol
   return Generic.PowerSeriesRing(R, prec, s; cached=cached, model=model)
end

function PowerSeriesRing(R::Ring, prec::Int,
   s::Vector{T}; cached=true, model=:capped_absolute) where
                                               T <: Union{Char, AbstractString}
   sym = [Symbol(v) for v in s]
   return Generic.PowerSeriesRing(R, prec, sym; cached=cached, model=model)
end
