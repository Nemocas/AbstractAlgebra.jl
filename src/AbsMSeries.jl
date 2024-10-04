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

@doc raw"""
    symbols(R::MSeriesRing)

Return a vector of symbols, one for each of the variables of the series ring
$R$.
"""
symbols(R::MSeriesRing) = R.sym

parent(a::MSeriesElem) = a.parent

base_ring_type(::Type{<:MSeriesRing{T}}) where T <: RingElement = parent_type(T)

function base_ring(R::MSeriesRing{T}) where T <: RingElement
    return base_ring(poly_ring(R))::parent_type(T)
end

@doc raw"""
    characteristic(a::MSeriesRing)

Return the characteristic of the base ring of the series `a`. If the
characteristic is not known, an exception is raised.
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
   R = parent(a)
   apoly = poly(a)

   poly_sum = Expr(:call, :+)
   n = nvars(parent(apoly))

   iter = zip(coefficients(apoly), exponent_vectors(apoly))
   citer = collect(iter)
   if R.weighted_prec != -1
      cv = sort!(citer; by=(tup->Base.sum(weights(R) .* tup[2])))
   else
      cv = reverse!(citer)
   end

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

   wp = parent(a).weighted_prec
   if wp == -1
      for i in nvars(parent(a)):-1:1
         push!(sum.args, Expr(:call, :O, Expr(:call, :^, x[i], a.prec[i])))
      end
   else
      push!(sum.args, Expr(:call, :O, :($wp)))
   end

   return sum
end

@enable_all_show_via_expressify MSeriesElem

function show(io::IO, mime::MIME"text/plain", p::MSeriesRing)
  @show_name(io, p)
  @show_special(io, mime, p)

  max_vars = 5 # largest number of variables to print
  n = nvars(p)
  print(io, "Multivariate power series ring")
  print(io, " in ", ItemQuantity(nvars(p), "variable"), " ")
  if n > max_vars
    join(io, symbols(p)[1:max_vars - 1], ", ")
    println(io, ", ..., ", symbols(p)[n])
  else
    join(io, symbols(p), ", ")
    println(io)
  end
  io = pretty(io)
  print(io, Indent(), "over ", Lowercase(), base_ring(p))
  print(io, Dedent())
end

function show(io::IO, p::MSeriesRing)
  @show_name(io, p)
  @show_special(io, p)
  if is_terse(io)
    print(io, "Multivariate power series ring")
  else
    io = pretty(io)
    print(io, "Multivariate power series ring in ", ItemQuantity(nvars(p), "variable"))
    print(terse(io), " over ", Lowercase(), base_ring(p))
  end
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::MSeriesRing, _, _) = elem_type(S)

function RandomExtensions.make(S::MSeriesRing,
                                      term_range::AbstractUnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, vs[1])
   else
      Make(S, term_range, make(R, vs...))
   end
end

function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make3{
                  <:RingElement, <:MSeriesRing, <:AbstractUnitRange{Int}}})
   S, term_range, v = sp[][1:end]
   f = S()
   g = gens(S)
   R = base_ring(S)
   if S.weighted_prec == -1
      prec = max_precision(S)
      for i = 1:rand(rng, term_range)
         term = S(1)
         for j = 1:length(g)
            term *= g[j]^rand(rng, 0:prec[j])
         end
         term *= rand(rng, v)
         f += term
      end
   else
      wt = weights(S)
      for i = 1:rand(rng, term_range)
         total = rand(0:S.weighted_prec)
         vv = Int[rand(0:total) for i = 1:length(g) - 1]
         vv = vcat(0, sort!(vv), total)
         w = Int[vv[i + 1] - vv[i] for i = 1:length(vv) - 1]
         ex = [Int(round(w[i]/wt[i])) for i in 1:length(w)]
         term = S(1)
         for j = 1:length(g)
            term *= g[j]^ex[j]
         end
         term *= rand(rng, v)
         f += term
      end
   end
   return f
end

function rand(rng::AbstractRNG, S::MSeriesRing,
                                             term_range::AbstractUnitRange{Int}, v...)
   rand(rng, make(S, term_range, v...))
end

@doc raw"""
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

@varnames_interface Generic.power_series_ring(R::Ring, prec::Int, s)
@varnames_interface Generic.power_series_ring(R::Ring, weights::Vector{Int}, prec::Int, s) macros=:no # use keyword `weights=...` instead
@varnames_interface Generic.power_series_ring(R::Ring, prec::Vector{Int}, s) n=:no macros=:no # `n` variant would clash with line above; macro would be the same as for `prec::Int`    
