###############################################################################
#
#   OrePolyRings.jl: skew-commutative polynomial ring R<D;σ,δ>
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

#coefficient_ring_type(T::Type{<:OrePolyRing}) = base_ring_type(T)
#coefficient_ring(R::OrePolyRing{T}) where T <: RingElement = base_ring(R)

#is_domain_type(::Type{OrePolyRingElem{T}}) where {T <: RingElement} = is_domain_type(T)
#is_exact_type(::Type{OrePolyRingElem{T}}) where {T <: RingElement} = is_exact_type(T)

number_of_variables(a::OrePolyRing) = 1

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::OrePolyRingElem, D=parent(a)|>var; context=nothing)
  sum = Expr(:call,:+)
  for i in order(a):-1:0
    c = coeff(a,i)
    if !iszero(c)
      xi = i < 1 ? 1 : i==1 ? D : Expr(:call,:^,D,i)
      if isone(c)
        push!(sum.args, Expr(:call,:*,xi))
      else
        push!(sum.args, Expr(:call,:*,expressify(c; context = context), xi))
      end
    end
  end

  return sum
end
@enable_all_show_via_expressify OrePolyRingElem

function show(io::IO, R::OrePolyRing)
  @show_name(io, R)
  @show_special(io, R)

  if is_terse(io)
    print(io, "Univariate Ore extension")
  else
    io = pretty(io)
    println(io, "Univariate Ore extension in ", gen(R))
    print(io,Indent(),"over ",Lowercase(),base_ring(R))
  end
end


###############################################################################
#
#   SkewDerivations
#
###############################################################################

function sigma_endomorphism(::SkewDerivation{D,S}) where {D,S} end
