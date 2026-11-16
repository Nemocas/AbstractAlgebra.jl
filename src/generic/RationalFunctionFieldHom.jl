################################################################################
#
#  Field access and printing
#
################################################################################

domain(f::RationalFunctionFieldHom) = f.domain

codomain(f::RationalFunctionFieldHom) = codomain(f.polynomial_map)

_coefficient_map(f::RationalFunctionFieldHom) = _coefficient_map(f.polynomial_map)

_image(f::RationalFunctionFieldHom) = _image(f.polynomial_map)

function Base.show(io::IO, f::RationalFunctionFieldHom)
  if is_terse(io)
    print(io, "Field homomorphism")
  else
    io = pretty(io)
    print(io, "Hom: ")
    print(terse(io), Lowercase(), domain(f), " -> ")
    print(terse(io), Lowercase(), codomain(f))
  end
end

function AbstractAlgebra.show_map_data(io::IO, f::RationalFunctionFieldHom)
  AbstractAlgebra.show_map_data(io, f.polynomial_map)
end

################################################################################
#
#  Types
#
################################################################################

function morphism_type(::Type{D}, ::Type{C}, ::Type{S}) where
    {T, U <: PolyRingElem{T}, D <: RationalFunctionField{T, U}, C <: Field, S}
  V = S === typeof(identity) || S <: IdentityMap ? Nothing : S
  P = PolyRingAnyMap{parent_type(U), C, V, elem_type(C)}
  return RationalFunctionFieldHom{D, C, P}
end

function morphism_type(::Type{D}, ::Type{C}) where
    {T, U <: PolyRingElem{T}, D <: RationalFunctionField{T, U}, C <: Field}
  return morphism_type(D, C, Nothing)
end

function morphism_type(F::RationalFunctionField{T, U}, E::Field, sigma) where
    {T, U <: PolyRingElem{T}}
  return morphism_type(typeof(F), typeof(E), typeof(sigma))
end

function morphism_type(F::RationalFunctionField{T, U}, E::Field) where
    {T, U <: PolyRingElem{T}}
  return morphism_type(typeof(F), typeof(E))
end

################################################################################
#
#  Construction
#
################################################################################

# There is no general transcendence test for arbitrary field implementations.
# Check the extension condition only when the available information permits it.
function _check_rational_function_field_hom(F, E, sigma, a)
  if is_known(is_finite, E)
    @req !is_finite(E) "A rational function field cannot embed into a finite field"
  end
  K = base_ring(F)
  @req sigma !== nothing || K !== E "The generator image is algebraic over the coefficient image"
  if E isa RationalFunctionField{<:Any, <:PolyRingElem}
    L = base_ring(E)
    n, d = _reduce_fraction(numerator(a, false), denominator(a, false))
    nonconstant = !is_constant(n) || !is_constant(d)
    if sigma === nothing && K === L
      @req nonconstant "The generator image must be transcendental over the coefficient field"
      return nothing
    elseif sigma isa Map && codomain(sigma) === L && nonconstant
      return nothing
    end
  end
  return nothing
end

@doc raw"""
    hom(F::Generic.RationalFunctionField, E::Field, [sigma,] a; check::Bool=true)

Return the homomorphism `F = K(t) -> E` sending `t` to `a` and acting
on `K` by the map `sigma` (canonical map if omitted).

The image `a` must be transcendental over `sigma(K)`.

Set `check=false` to skip validation; with `check=true`, only available checks
are performed.
"""
function hom(F::RationalFunctionField{T, U}, E::Field, sigma, a; check::Bool=true) where
    {T, U <: PolyRingElem{T}}
  K = base_ring(F)
  if check
    if sigma isa Map
      @req domain(sigma) === K "The coefficient map must have domain equal to the coefficient field"
      @req codomain(sigma) isa Field "The coefficient map must have a field codomain"
    end
    if is_known(characteristic, F) && is_known(characteristic, E)
      @req characteristic(F) == characteristic(E) "The fields must have the same characteristic"
    end
    if a isa elem_type(E)
      @req parent(a) === E "The generator image must belong to the codomain"
    end
  end
  if sigma === identity || sigma isa IdentityMap
    sigma = nothing
  end
  img = E(a)
  if check
    @req parent(img) === E "The generator image must belong to the codomain"
    _check_rational_function_field_hom(F, E, sigma, img)
  end
  R = base_ring(underlying_fraction_field(F))
  p = hom(R, E, sigma, img)
  return RationalFunctionFieldHom{typeof(F), typeof(E), typeof(p)}(F, p)
end

function hom(F::RationalFunctionField{T, U}, E::Field, a; check::Bool=true) where
    {T, U <: PolyRingElem{T}}
  return hom(F, E, nothing, a; check=check)
end

################################################################################
#
#  Evaluation and composition
#
################################################################################

function image(f::RationalFunctionFieldHom, x)
  F = domain(f)
  if x isa elem_type(F)
    @req parent(x) === F "Element not in the domain"
  else
    x = F(x)
  end
  E = codomain(f)
  p = f.polynomial_map
  n = E(p(numerator(x, false)))
  d = E(p(denominator(x, false)))
  @req !iszero(d) "The map is illegal because the denominator maps to zero"
  return divexact(n, d)::elem_type(E)
end

(f::RationalFunctionFieldHom)(x) = image(f, x)

function compose(f::RationalFunctionFieldHom{D, C}, g::RationalFunctionFieldHom{C, E}) where {D, C, E}
  check_composable(f, g)
  F = domain(f)
  a = g(_image(f))
  sigma = _coefficient_map(f)
  tau = _coefficient_map(g)
  L = base_ring(domain(g))
  if sigma === nothing && base_ring(F) === L
    return hom(F, codomain(g), tau, a; check=false)
  elseif sigma isa Map && codomain(sigma) === L
    if tau === nothing
      return hom(F, codomain(g), sigma, a; check=false)
    elseif tau isa Map
      return hom(F, codomain(g), compose(sigma, tau), a; check=false)
    end
  end
  # Coefficient images can involve the intermediate generator, so in general
  # the second map must act on their images in the full intermediate field.
  sigma = map_from_func(base_ring(F), codomain(g), c -> g(f(F(c))))
  return hom(F, codomain(g), sigma, a; check=false)
end

################################################################################
#
#  Möbius inversion
#
################################################################################

function Base.inv(f::RationalFunctionFieldHom)
  F = domain(f)
  E = codomain(f)
  @req F === E "inv only implemented for endomorphisms"
  n, d = _reduce_fraction(numerator(_image(f), false), denominator(_image(f), false))
  # [K(t):K(n/d)] = max(deg(n),deg(d)) for coprime n,d; see
  # Stichtenoth, Algebraic Function Fields and Codes, 2nd ed. (Springer, 2009),
  # Theorem 1.4.11 and Example 1.4.13, pp. 19–21.
  @req max(degree(n), degree(d)) == 1 "The generator image must be a Möbius transformation"
  a, b = coeff(n, 1), coeff(n, 0)
  c, e = coeff(d, 1), coeff(d, 0)
  sigma = _coefficient_map(f)
  if sigma === nothing
    sigma_inv = nothing
  else
    @req sigma isa Map && codomain(sigma) === base_ring(E) "Möbius inversion requires an invertible map between the coefficient fields"
    if sigma isa AbstractAlgebra.MapFromFunc
      @req isdefined(sigma, :g) "An inverse coefficient map is required"
    end
    sigma_inv = inv(sigma)
    @req domain(sigma_inv) === base_ring(E) && codomain(sigma_inv) === base_ring(F) "The inverse coefficient map has incompatible domain or codomain"
    a, b, c, e = sigma_inv(a), sigma_inv(b), sigma_inv(c), sigma_inv(e)
  end
  t = gen(F)
  img = divexact(F(e)*t - F(b), F(a) - F(c)*t)
  return hom(E, F, sigma_inv, img; check=false)
end

function preimage(f::RationalFunctionFieldHom, y)
  @req parent(y) === codomain(f) "Element not in the codomain"
  return inv(f)(y)
end
