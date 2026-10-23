################################################################################
#
#  Field access
#
################################################################################
domain(f::PolyFracFieldAnyMap{D,C,V}) where {D,C,V} = f.domain
codomain(f::PolyFracFieldAnyMap{D,C,V}) where {D,C,V} = f.morphism.codomain
underlying_morphism(f) = f.morphism

################################################################################
#
#  String I/O
#
################################################################################
function Base.show(io::IO, f::PolyFracFieldAnyMap)
  io = pretty(io)
  if is_terse(io)
    print(io, "Ring homomorphism")
  else
    print(io, "Hom: ")
    print(terse(io), Lowercase(), domain(f), " -> ")
    print(terse(io), Lowercase(), codomain(f))
  end
end

function AbstractAlgebra.show_map_data(io::IO, f::PolyFracFieldAnyMap)
  println(io)
  println(io, "defined by", Indent())
  R = domain(f)
  g = gen(R)
  print(io, g, " -> ", f(g))
  print(io, Dedent())
end

################################################################################
#
#  Constructor
#
################################################################################
function hom(R::D, S::NCRing, image) where {D<:FracField{<:PolyRingElem}}
  r = base_ring(R)
  return PolyFracFieldAnyMap(R,hom(r,S,image))
end

function hom(R::D, S::NCRing, image) where {D<:RationalFunctionField{<:RingElement,<:PolyRingElem}}
  r = R() |> numerator |> parent
  return PolyFracFieldAnyMap(R,hom(r,S,image))
end

################################################################################
#
#  Evaluation functions
#
################################################################################
function (f::PolyFracFieldAnyMap{D,C,V})(a) where {D,C,V<:Map{<:PolyRing,C}}
  ϕ = underlying_morphism(f)
  return ϕ.(numerator(a))/ϕ.(denominator(a))
end
