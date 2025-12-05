################################################################################
#
#  Field access
#
################################################################################

domain(f::PolyRingAnyMap) = f.domain
codomain(f::PolyRingAnyMap) = f.codomain

# Not sure if we want to expose the following function to the user.
# It might be `nothing`. We could return `identity` in the `nothing` case.
_coefficient_map(f::PolyRingAnyMap) = f.coeff_map

_image(f::PolyRingAnyMap) = f.img_gen

################################################################################
#
#  String I/O
#
################################################################################

function AbstractAlgebra.show_map_data(io::IO, f::PolyRingAnyMap)
  println(io)
  println(io, "defined by", Indent())
  R = domain(f)
  g = gen(R)
  print(io, g, " -> ", f(g))
  print(io, Dedent())
  # the last print statement must not add a new line
  phi = _coefficient_map(f)
  if !(phi isa Nothing)
    println(io)
    println(io, "with map on coefficients")
    print(io, Indent(), phi, Dedent())
  end
end

function Base.show(io::IO, f::PolyRingAnyMap)
  io = pretty(io)
  if is_terse(io)
    print(io, "Ring homomorphism")
  else
    print(io, "Hom: ")
    print(terse(io), Lowercase(), domain(f), " -> ")
    print(terse(io), Lowercase(), codomain(f))
  end
end

################################################################################
#
#  Helper
#
################################################################################

# # Since we want to allow specifying images in a "subring", we need to coerce
# if necessary. For example, hom(Qx, Qx, [1, 1]), should work, although
# 1 is not an element of the codomain.
function _coerce(S, img_gen)
  if typeof(img_gen) === elem_type(S)
    return img_gen::elem_type(S)
  else
    _img_gen = S(img_gen)
    typeof(_img_gen) === elem_type(S) || error("Elements cannot be coerced into the codomain")
    return _img_gen::elem_type(S)
  end
end

# When evaluating the map F at a polynomial f, we first construct the polynomial
# map_coefficients(_coefficient_map(F), f), which is a polynomial over
# codomain(_coefficient_map(F)).

function temp_ring(f::PolyRingAnyMap{<:Any, <: Any, <: Map})
  if isdefined(f, :temp_ring)
    return f.temp_ring::poly_ring_type(codomain(_coefficient_map(f)))
  end

  S, = polynomial_ring(codomain(_coefficient_map(f)), cached = false)
  f.temp_ring = S
  return S
end

# If the _coefficient_map is e.g. a julia function, there is not too much we can
# do, so we do the defensive thing
function temp_ring(f::PolyRingAnyMap{<:Any, <: Any})
  return nothing
end
################################################################################
#
#  Composition
#
################################################################################

# This is getting difficult, because Map{C, D} does not yield information
# on the type of the domain, codomain

# First consider the case where both coefficient maps are maps in the Map
# sense
function compose(F::PolyRingAnyMap{D, C, S}, G::PolyRingAnyMap{C, E, U}) where {D, C, E, S <: Map, U <: Map}
  codomain(F) === domain(G) || error("Incompatible (co)domain in composition")
  f = _coefficient_map(F)
  g = _coefficient_map(G)
  if typeof(codomain(f)) === typeof(domain(g))
    newcoeffmap = compose(f, g)
    return hom(domain(F), codomain(G), newcoeffmap, G.(_image(F)))
  else
    return Generic.CompositeMap(F, G)
  end
end

# No coefficient maps in both maps
function compose(F::PolyRingAnyMap{D, C, Nothing}, G::PolyRingAnyMap{C, E, Nothing}) where {D, C, E}
  codomain(F) === domain(G) || error("Incompatible (co)domain in composition")
  return hom(domain(F), codomain(G), G.(_image(F)))
end

# Julia functions in both maps
function compose(F::PolyRingAnyMap{D, C, <: Function}, G::PolyRingAnyMap{C, E, <: Function}) where {D, C, E}
  codomain(F) === domain(G) || error("Incompatible (co)domain in composition")
  b = _coefficient_map(F)(one(coefficient_ring(domain(F))))
  if parent(b) === domain(G)
    return hom(domain(F), codomain(G), x -> G(_coefficient_map(F)(x)), G.(_image(F)))
  elseif parent(b) === coefficient_ring(domain(G))
    return hom(domain(F), codomain(G), x -> _coefficient_map(G)(_coefficient_map(F)(x)), G.(_image(F)))
  else
    error("coefficient map is not admissible")
  end
end

# Now compose with arbitrary maps
# For now we only properly do the Map version of the coefficient map has codomain C
# Otherwise we bail out and return a generic CompositeMap
function compose(F::PolyRingAnyMap{D, C, <: Map, <: Any}, G::S) where {D, C, S <: Map{C, <: Any}}
  codomain(F) === domain(G) || error("Incompatible (co)domain in composition")
  f = _coefficient_map(F)
  if typeof(codomain(f)) === C
    newcoeffmap = compose(f, G)
    return hom(domain(F), codomain(G), newcoeffmap, G.(_image(F)))
  else
    return Generic.CompositeMap(F, G)
  end
end

# Special case for composing with the identity
function compose(F::PolyRingAnyMap{D, C, <: Map, <: Any}, G::S) where {D, C, S <: Generic.IdentityMap{C}}
  codomain(F) === domain(G) || error("Incompatible (co)domain in composition")
  f = _coefficient_map(F)
  if typeof(codomain(f)) === C
    newcoeffmap = compose(f, G)
    return hom(domain(F), codomain(G), newcoeffmap, G.(_image(F)))
  else
    return Generic.CompositeMap(F, G)
  end
end

################################################################################
#
#  Types computers
#
################################################################################

function morphism_type(::Type{D}, ::Type{C}) where {D <: AbstractAlgebra.PolyRing, C <: NCRing}
  return PolyRingAnyMap{D, C, Nothing, elem_type(C)}
end

morphism_type(::D, ::C) where {D <: AbstractAlgebra.PolyRing, C <: NCRing} = morphism_type(D, C)

function morphism_type(::Type{D}, ::Type{C}, f::Type{F}) where {D <: AbstractAlgebra.PolyRing, C <: NCRing, F}
  return PolyRingAnyMap{D, C, F, elem_type(C)}
end

morphism_type(::D, ::C, ::F) where {D <: AbstractAlgebra.PolyRing, C <: NCRing, F} = morphism_type(D, C, F)

################################################################################
#
#  Constructor
#
################################################################################

@doc raw"""
    hom(R::AbstractAlgebra.PolyRing, S::NCRing, [coeff_map,] image)
    
Given a homomorphism `coeff_map` from `C` to `S`, where `C` is the 
coefficient ring of `R`, and given an element `image` of `S`, return the
homomorphism from `R` to `S` whose restriction 
to `C` is `coeff_map`, and which sends the generator of `R` to `image`.
 
If no coefficient map is entered, invoke a canonical homomorphism of `C`
to `S`, if such a homomorphism exists, and throw an error, otherwise.

# Examples
```jldoctest
julia> Zx, x = ZZ[:x];

julia> F = hom(Zx, Zx, x + 1);

julia> F(x^2)
x^2 + 2*x + 1

julia> Fp = GF(3); Fpy, y = Fp[:y];

julia> G = hom(Zx, Fpy, c -> Fp(c), y^3);

julia> G(5*x + 1)
2*y^3 + 1
```
"""
function hom(R::AbstractAlgebra.PolyRing, S::NCRing, coeff_map, image)
  # Now coerce into S or throw an error if not possible
  img = _coerce(S, image)
  return PolyRingAnyMap(R, S, coeff_map, img)
end

function hom(R::AbstractAlgebra.PolyRing, S::NCRing, image)
  # Now coerce into S or throw an error if not possible
  img = _coerce(S, image)
  return PolyRingAnyMap(R, S, nothing, img)
end

################################################################################
#
#  Evaluation functions
#
################################################################################

# if the coefficient map is nothing
function _evaluate_help(F::PolyRingAnyMap{<: AbstractAlgebra.PolyRing, <: Any, Nothing}, g)
  return g(F.img_gen)
end

# general case
function _evaluate_help(F::PolyRingAnyMap{<: AbstractAlgebra.PolyRing, <: Any, <: Any}, g)
  @assert !(_coefficient_map(F) isa Nothing)
  S = temp_ring(F)
  if S !== nothing
    return (map_coefficients(_coefficient_map(F), g, parent = S))(F.img_gen)
  else
    return (map_coefficients(_coefficient_map(F), g))(F.img_gen)
  end
end

function (F::PolyRingAnyMap)(g)
  if g isa elem_type(domain(F))
    _evaluate_help(F, g)
  else 
    gg = domain(F)(g)
    @assert parent(gg) === domain(F)
    return F(gg)
  end
end
