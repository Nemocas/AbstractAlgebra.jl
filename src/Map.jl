################################################################################
#
#  Map.jl : Maps
#
################################################################################

################################################################################
#
#  All maps
#
################################################################################

Map(::Type{T}) where T <: Map = supertype(T)
Map(::Type{S}) where S <: SetMap = Map{D, C, <:S, T} where {D, C, T}

function domain end
function codomain end
function image_fn end

function coimage(h::Map)
  return quo(domain(h), kernel(h)[1])
end

function check_composable(a::Map, b::Map)
   codomain(a) !== domain(b) && error("Incompatible maps")
end

function Base.show(io::IO, mime::MIME"text/plain", M::Map)
   # the "header" printed for most map types is identical: a descriptive name
   # for the map followed by "from DOMAIN" and "to CODOMAIN" lines. We
   # implement this generically here. If desired, map types can provide a
   # custom show_map_data method to add additional data after this initial
   # "header".
   #
   # The output will be
   #
   # show_map_head(io, M)
   #   from DOMAIN
   #   to CODOMAIN
   # show_map_data(io, M)
   #
   show_map_head(io, M)
   println(io)
   io = pretty(io)
   println(io, Indent(), "from ", Lowercase(), domain(M))
   print(io, "to ", Lowercase(), codomain(M), Dedent())
   show_map_data(io, M)
end

# for backwards compatibility, but all maps using `@show_name` should overload it
show_map_head(io::IO, M::Map) = print(terse(io), M)

function show_map_data(io::IO, M::Map)
   # no extra data by default; Map subtypes may overload this
end

function Base.show(io::IO, M::Map)
   if is_terse(io)
      print(io, "Map")
   else
      io = pretty(io)
      io = terse(io)
      print(io, "Map: ")
      print(io, Lowercase(), domain(M), " -> ")
      print(io, Lowercase(), codomain(M))
   end
end

Base.broadcastable(M::Map) = Ref(M)

Base.:\(f::Map, x) = preimage(f, x)

###############################################################################
#
#   CompositeMap
#
###############################################################################

@doc raw"""
    compose(f::Map, g::Map)

Compose the two maps $f$ and $g$, i.e. return the map $h$ such that $h(x) = g(f(x))$.

# Examples
```jldoctest
julia> f = map_from_func(ZZ, ZZ, x -> x + 1);

julia> g = map_from_func(ZZ, QQ, x -> QQ(x));

julia> h = compose(f, g)
Composite map
  from integers
  to rationals
which is the composite of
  Map: integers -> integers
  Map: integers -> rationals
```
"""
function compose(f::Map, g::Map)
   check_composable(f, g)
   return Generic.CompositeMap(f, g)
end

*(f::Map, g::Map) = compose(f, g)

################################################################################
#
#  IdentityMap
#
################################################################################

@doc raw"""
    identity_map(R::D) where D <: AbstractAlgebra.Set

Return an identity map on the domain $R$.

# Examples
```jldoctest
julia> R, t = ZZ[:t]
(Univariate polynomial ring in t over integers, t)

julia> f = identity_map(R)
Identity map
  of univariate polynomial ring in t over integers

julia> f(t)
t
```
"""
identity_map(R::D) where D <: Set = Generic.IdentityMap{D}(R)

################################################################################
#
#  FunctionalMap
#
################################################################################

@doc raw"""
    map_from_func(D, C, image_fn, [inverse_fn])

Creates the map `D -> C, x -> image_fn(x)` of type `MapFromFunc` given the callable
object `image_fn`. If `inverse_fn` is provided, it is assumed to satisfy
`image_fn(inverse_fn(x)) = x` and will be used as the preimage function.

Construct the MapFromFunc with domain and codomain given by the parent objects
$R$ and $S$ corresponding to the Julia function $f$.

# Examples
```jldoctest
julia> f = map_from_func(ZZ, ZZ, x -> x + 1)
Map defined by a Julia function
  from integers
  to integers

julia> f(ZZ(2))
3
```
"""
map_from_func(D, C, image_fn) = MapFromFunc(D, C, image_fn)
map_from_func(D, C, image_fn, inverse_fn) = MapFromFunc(D, C, image_fn, inverse_fn)


################################################################################
#
#  HeckeMap
#
################################################################################


# Hecke maps store attributes in the header object
_get_attributes(G::Map{<:Any, <:Any, HeckeMap, <:Any}) = _get_attributes(G.header)
_get_attributes!(G::Map{<:Any, <:Any, HeckeMap, <:Any}) = _get_attributes!(G.header)
_is_attribute_storing_type(::Type{<:Map(HeckeMap)}) = true

(f::Map(HeckeMap))(x) = image(f, x)

function domain(M::Map(HeckeMap))
  return M.header.domain
end

function codomain(M::Map(HeckeMap))
  return M.header.codomain
end

function image_function(f::Map(HeckeMap))
  if isdefined(f.header, :image)
    return f.header.image
  else
    return x -> image(f, x)
  end
end

function preimage_function(f::Map(HeckeMap))
  if isdefined(f.header, :preimage)
    return f.header.preimage
  else
    return x -> preimage(f, x)
  end
end

image_fn(f::Map(HeckeMap)) = image_function(f)
preimage_fn(f::Map(HeckeMap)) = preimage_function(f)


################################################################################
#
#  Comparison of objects as maps
#
#  Rationale: Often objects are implicitly used as maps. For instance, a `Ring`
#  `R` can be used to type-cast objects `x` into `R`. When objects are
#  implicitly used as maps, one may wish to compare them as such. The function
#  `==` can not be overloaded for this since it would be expected to compare
#  two objects as what they are originally (i.e. as `Ring`s in the above
#  example). Therefore, we introduce another function here for which methods
#  can be added successively to allow for reasonable comparisons.
#
#  This is, for example, used in recursive setups like the following.
#  Suppose `R = ℤ[x₁, x₂]` and `P = ℚ[y₁, y₂]`, together with maps
#
#    φ, ψ : R → P, xᵢ ↦ yᵢ
#
#  Both maps will need to implicitly cast the coefficients in `R` to
#  coefficients in `P`. Let us now assume that `φ` is doing this implicitly,
#  while `ψ` has `QQ` stored internally as a `coefficient_map`.
#    Comparing those two maps (i.e. implementing `==`) will now compare the
#  images of the variables (which agree in this case) and then move on to
#  comparing the `coefficient_map`s. In this case, this would compare
#  `nothing` to `QQ`. Clearly, `nothing == QQ` should return `false` by any
#  means. With the function `is_equal_as_morphism` introduced below, we
#  could implement a reasonable method for this case.
################################################################################

function is_equal_as_morphism(a::Any, b::Any)
  a === b && return true
  error("no method implemented to compare $a and $b as morphisms beyond `===`")
end
