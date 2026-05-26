
function trivial_derivation(R::D,σ::S) where {D<:Univariateish,S}
  return PolySkewDerivation{D,S}(R,σ,zero(R))
end

function derivation(R::D,σ::S,c::T) where {D<:Univariateish,S<:Map{D,D},T}
  @req elem_type(R) == T "incompatible coefficient type"
  return PolySkewDerivation{D,S}(R,σ,c)
end
derivation(R::D) where {D<:Univariateish} = derivation(R,identity_map(R),one(R))
derivation(R::D,σ::S) where {D<:Univariateish,S<:Map{D,D}} = derivation(R,σ,one(R))

function show(io::IO, d::PolySkewDerivation)
  io = pretty(io)
  if is_terse(io)
    print(io, "Skew-derivation")
  else
    print(io, "σ-Der: ")
    print(terse(io), Lowercase(), domain(d), " -> ")
    print(terse(io), Lowercase(), codomain(d))
  end
end

function show(io::IO, d::PolySkewDerivation{D,<:Map(IdentityMap)}) where D
  io = pretty(io)
  if is_terse(io)
    print(io, "Derivation")
  else
    print(io, "Der: ")
    print(terse(io), Lowercase(), domain(d), " -> ")
    print(terse(io), Lowercase(), codomain(d))
  end
end

function show_map_data(io::IO, d::PolySkewDerivation)
  println(io)
  if coefficient(d) |> iszero
    print(io, "which is identically zero ")
    println(io, "and with commutation rule", Indent())
  elseif !(coefficient(d) |> isone)
    println(io, "with coefficient ", coefficient(d))
    println(io, "and commutation rule", Indent())
  else
    println(io, "with commutation rule", Indent())
  end
  print(io, sigma_endomorphism(d))
  show_map_data(io,sigma_endomorphism(d))
end

function show_map_data(io::IO, d::PolySkewDerivation{D,<:Map(IdentityMap)}) where D
  if !(coefficient(d) |> isone)
    println(io)
    print(io, "with coefficient ", coefficient(d))
  end
end

domain(d::PolySkewDerivation) = d.domain
codomain(d::PolySkewDerivation) = d.domain
coefficient(d::PolySkewDerivation) = first(d.intermediate_cache)

sigma_endomorphism(δ::PolySkewDerivation{D,S}) where {D,S} = δ.σ

function (d::PolySkewDerivation{D,S})(a) where {D,S<:Map(IdentityMap)}
  c = coefficient(d)
  iszero(c) && return domain(d)()
  return c*derivative(a)
end

function (d::PolySkewDerivation{D,S})(a::T) where {D,S<:Map{D,D,<:Map,<:Any},T<:PolyRingElem}
  c = coefficient(d)
  iszero(c) && return domain(d)()

  x = parent(a) |> gen
  σ = sigma_endomorphism(d)
  cached_degree = length(d.intermediate_cache)+1
  for i in cached_degree:degree(a)
    res = last(d.intermediate_cache)*x + σ(x^(i-1))*c
    push!(d.intermediate_cache, res)
  end

  return sum(σ.(Iterators.drop(coefficients(a),1)) .* d.intermediate_cache; init=parent(a)())
end

function (d::PolySkewDerivation{D,S})(a::T) where {
  D<:FracField{<:PolyRingElem},
  S<:Map{D,D,<:Map,<:Any},
  T<:FracElem{<:PolyRingElem}}
  
  c = coefficient(d)
  iszero(c) && return domain(d)()

  R = codomain(d)
  σ = sigma_endomorphism(d)
  p = numerator(a)
  q = denominator(a)
  return R((d(p) - σ(q)*d(p))//(σ(q)*q))::elem_type(R)
end

function (d::PolySkewDerivation{D,S})(a::T) where {
  D<:RationalFunctionField{<:RingElement,<:PolyRingElem},
  S<:Map{D,D,<:Map,<:Any},
  T<:RationalFunctionFieldElem{<:RingElement,<:PolyRingElem}}

  c = coefficient(d)
  iszero(c) && return domain(d)()

  R = codomain(d)
  σ = sigma_endomorphism(d)
  p = numerator(a)
  q = denominator(a)
  return R((d(p) - σ(q)*d(p))//(σ(q)*q))::elem_type(codomain(d))
end
