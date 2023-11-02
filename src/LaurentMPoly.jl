###############################################################################
#
#   LaurentMPoly.jl: Multivariate Laurent polynomials
#
###############################################################################

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::LaurentMPolyRingElem, x = symbols(parent(a)); context = nothing)
    sum = Expr(:call, :+)
    n = nvars(parent(a))
    for (c, v) in zip(coefficients(a), exponent_vectors(a))
        prod = Expr(:call, :*)
        if !isone(c)
            push!(prod.args, expressify(c, context = context))
        end
        for i in 1:n
            if isone(v[i])
                push!(prod.args, x[i])
            elseif !iszero(v[i])
                push!(prod.args, Expr(:call, :^, x[i], v[i]))
            end
        end
        push!(sum.args, prod)
    end
    return sum
end

@enable_all_show_via_expressify LaurentMPolyRingElem

function show(io::IO, ::MIME"text/plain", p::LaurentMPolyRing)
  max_vars = 5 # largest number of variables to print
  n = nvars(p)
  print(io, "Multivariate Laurent polynomial ring")
  print(io, "in ", ItemQuantity(nvars(p), "variable"), " ")
  if n > max_vars
    join(io, symbols(p)[1:max_vars - 1], ", ")
    println(io, "..., ", symbols(p)[n])
  else
    join(io, symbols(p), ", ")
    println(io)
  end
  io = pretty(io)
  print(io, Indent(), "over ", Lowercase(), base_ring(p))
  print(io, Dedent())
end

function show(io::IO, p::LaurentMPolyRing)
  if get(io, :supercompact, false)
    # no nested printing
    print(io, "Multivariate Laurent polynomial ring")
  else
    # nested printing allowed, preferably supercompact
    io = pretty(io)
    print(io, "Multivariate Laurent polynomial ring in ", ItemQuantity(nvars(p), "variable"))
    print(IOContext(io, :supercompact => true), " over ", Lowercase(), base_ring(p))
  end
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function gens(R::LaurentMPolyRing)
    return [gen(R, i) for i in 1:nvars(R)]
end

coefficient_ring(a::LaurentMPolyRingElem) = coefficient_ring(parent(a))
base_ring(a::LaurentMPolyRingElem) = base_ring(parent(a))

###############################################################################
#
#   Derivative
#
###############################################################################

function derivative(a::LaurentMPolyRingElem{T}, x::LaurentMPolyRingElem{T}) where T <: RingElement
   check_parent(a, x)
   return derivative(a, var_index(x))
end

###############################################################################
#
#   Random elements
#
###############################################################################

function RandomExtensions.maketype(S::LaurentMPolyRing, _, _, _)
    return elem_type(S)
end

function RandomExtensions.make(S::LaurentMPolyRing,
                               term_range::AbstractUnitRange{Int},
                               exp_bound::AbstractUnitRange{Int},
                               vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, exp_bound, vs[1])
   else
      Make(S, term_range, exp_bound, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make4{<:RingElement,
                                         <:LaurentMPolyRing,
                                         <:AbstractUnitRange{Int},
                                         <:AbstractUnitRange{Int}}})
   S, term_range, exp_bound, v = sp[][1:end]
   f = zero(S)
   g = gens(S)
   R = base_ring(S)
   for i = 1:rand(rng, term_range)
      term = one(S)
      for j = 1:length(g)
         term *= g[j]^rand(rng, exp_bound)
      end
      term *= rand(rng, v)
      f += term
   end
   return f
end

function rand(rng::AbstractRNG, S::LaurentMPolyRing,
              term_range::AbstractUnitRange{Int}, exp_bound::AbstractUnitRange{Int}, v...)
   rand(rng, make(S, term_range, exp_bound, v...))
end

function rand(S::LaurentMPolyRing, term_range, exp_bound, v...)
   rand(GLOBAL_RNG, S, term_range, exp_bound, v...)
end

###############################################################################
#
#   laurent_polynomial_ring constructor
#
###############################################################################

@doc raw"""
    laurent_polynomial_ring(R::Ring, s::Vector{T}; cached::Bool = true) where T <: VarName

Given a base ring `R` and an array of strings `s` specifying how the
generators (variables) should be printed, return a tuple `T, (x1, x2, ...)`
representing the new ring $T = R[x1, 1/x1, x2, 1/x2, ...]$ and the generators
$x1, x2, ...$ of the  ring. By default the parent object `T` will
depend only on `R` and `x1, x2, ...` and will be cached. Setting the optional
argument `cached` to `false` will prevent the parent object `T` from being
cached.
"""
function laurent_polynomial_ring(R::Ring, s::AbstractVector{<:VarName}; cached::Bool = true)
   return Generic.laurent_polynomial_ring(R, [Symbol(v) for v in s], cached=cached)
end

function laurent_polynomial_ring(R::Ring, s::Vector{Symbol}; cached::Bool = true)
   return Generic.laurent_polynomial_ring(R, s; cached=cached)
end

function laurent_polynomial_ring(R::Ring, n::Int, s::VarName=:x; cached::Bool = false)
   return Generic.laurent_polynomial_ring(R, [Symbol(s, i) for i=1:n]; cached=cached)
end
