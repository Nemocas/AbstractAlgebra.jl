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

function check_composable(a::Map, b::Map)
   codomain(a) != domain(b) && error("Incompatible maps")
end

function Base.show(io::IO, ::MIME"text/plain", M::AbstractAlgebra.Map)
   # the "header" is identical to the supercompact output; this
   # allows other map types to reuse this method
   println(IOContext(io, :supercompact => true), M)
   io = pretty(io)
   println(io, Indent(), "from ", Lowercase(), domain(M))
   print(io, "to ", Lowercase(), codomain(M), Dedent())
   show_map_data(io, M)
end

function show_map_data(io::IO, M::AbstractAlgebra.Map)
   # no extra data by default; specific Map subtypes can overload this
end

function Base.show(io::IO, M::AbstractAlgebra.Map)
   if get(io, :supercompact, false)
      # no nested printing
      print(io, "Map")
   else
      # nested printing allowed, preferably supercompact
      io = pretty(io)
      print(io, "Map: ")
      print(IOContext(io, :supercompact => true), Lowercase(), domain(M), " -> ")
      print(IOContext(io, :supercompact => true), Lowercase(), codomain(M))
   end
end

###############################################################################
#
#   CompositeMap
#
###############################################################################

@doc raw"""
    compose(f::Map, g::Map)

Compose the two maps $f$ and $g$, i.e. return the map $h$ such that $h(x) = g(f(x))$.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> f = map_from_func(x -> x + 1, ZZ, ZZ);

julia> g = map_from_func(x -> QQ(x), ZZ, QQ);

julia> h = compose(f, g)
Functional composite map
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
```jldoctest; setup = :(using AbstractAlgebra)
julia> R, t = ZZ[:t]
(Univariate polynomial ring in t over integers, t)

julia> f = identity_map(R)
Identity map
  of univariate polynomial ring in t over integers

julia> f(t)
t
```
"""
identity_map(R::D) where D <: AbstractAlgebra.Set = Generic.IdentityMap{D}(R)

################################################################################
#
#  FunctionalMap
#
################################################################################

@doc raw"""
    map_from_func(image_fn::Function, domain, codomain)

Construct the generic functional map with domain and codomain given by the parent objects
$R$ and $S$ corresponding to the Julia function $f$.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> f = map_from_func(x -> x + 1, ZZ, ZZ)
Map defined by a Julia function
  from integers
  to integers

julia> f(ZZ(2))
3
```
"""
function map_from_func(image_fn::Function, domain, codomain)
   return Generic.FunctionalMap(domain, codomain, image_fn)
end

################################################################################
#
#  FunctionalCompositeMap
#
################################################################################

function compose(f::Map(FunctionalMap){D, U}, g::Map(FunctionalMap){U, C}) where {D, U, C}
   check_composable(f, g)
   return Generic.FunctionalCompositeMap(f, g)
end

