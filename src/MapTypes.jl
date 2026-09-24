###############################################################################
#
#  MapTypes.jl : Maps
#
###############################################################################

mutable struct MapCache{D, C, De, Ce}
  lim::Int

  im::Dict{De, Ce}
  imStat::Dict{De, Int}

  pr::Dict{Ce, De}
  prStat::Dict{Ce, Int}

  old_im
  old_pr

  function MapCache(dom::D, cod::C, ::Type{De}, ::Type{Ce}, lim::Int=100) where {D, C, De, Ce}
    r = new{D, C, De, Ce}()
    r.lim = lim
    r.im = Dict{De, Ce}()
    r.imStat = Dict{De, Int}()
    r.pr = Dict{Ce, De}()
    r.prStat = Dict{Ce, Int}()
    return r
  end
end

@attributes mutable struct MapHeader{D, C}
  domain::D
  codomain::C
  image           # object used to compute images
  preimage        # object used to compute preimages
  cache::MapCache

  function MapHeader{D, C}() where {D, C}
    z = new{D, C}()
    return z
  end

  function MapHeader{D, C}(domain::D, codomain::C) where {D, C}
    return new{D,C}(domain, codomain)
  end

  function MapHeader{D, C}(domain::D, codomain::C, image) where {D, C}
    return new{D, C}(domain, codomain, image)
  end

  function MapHeader{D, C}(domain::D, codomain::C, image, preimage) where {D, C}
    return new{D, C}(domain, codomain, image, preimage)
  end
end

function MapHeader(domain::D, codomain::C) where {D, C}
  return MapHeader{D, C}(domain, codomain)
end

function MapHeader(domain::D, codomain::C, image) where {D, C}
  return MapHeader{D, C}(domain, codomain, image)
end

function MapHeader(domain::D, codomain::C, image, preimage) where {D, C}
  return MapHeader{D, C}(domain, codomain, image, preimage)
end

# preimage_function(f) = a -> preimage(f, a)
# image_function(f) = a -> image(f, a)


###########################################################
# To turn a Function (method) into a map.
###########################################################

"""
    MapFromFunc(D, C, f, [g])

Creates the map `D -> C, x -> f(x)` given the callable
object `f`. If `g` is provided, it is assumed to satisfy
`f(g(x)) = x` and will be used as the preimage function.

# Example

```jldoctest
julia> F = GF(2);

julia> f = MapFromFunc(QQ, F, x -> F(numerator(x)) * inv(F(denominator(x))))
Map defined by a Julia function
  from rationals
  to finite field F_2

julia> f(QQ(1//3))
1

julia> println(f)
Map: rationals -> F

julia> f = MapFromFunc(QQ, F, x -> F(numerator(x)) * inv(F(denominator(x))), y -> QQ(lift(y)),)
Map defined by a Julia function with inverse
  from rationals
  to finite field F_2

julia> preimage(f, F(1))
1//1

julia> println(f)
Map: rationals -> F

```
"""
mutable struct MapFromFunc{R, T} <: Map{R, T, MapWithHeader, MapFromFunc}
  header::MapHeader{R, T}
  f
  g

  function MapFromFunc{R, T}(D::R, C::T, f) where {R, T}
    n = new{R, T}()
    n.header = MapHeader(D, C, f)
    n.f = f
    return n
  end

  function MapFromFunc{R, T}(D::R, C::T, f, g) where {R, T}
    n = new{R, T}()
    n.header = MapHeader(D, C, f, g)
    n.f = f
    n.g = g
    return n
  end
end

function image(f::MapFromFunc, x)
  @req parent(x) === domain(f) "Element not in the domain"
  y = f.f(x)
  @req parent(y) === codomain(f) "Image not in the codomain"
  return y::elem_type(codomain(f))
end

function preimage(f::MapFromFunc, y)
  @req parent(y) === codomain(f) "Element not in the codomain"
  x = f.g(y)::elem_type(domain(f))
  @req parent(x) === domain(f) "Preimage not in the domain"
  return x
end

function Base.show(io::IO, ::MIME"text/plain", M::MapFromFunc)
  @show_name(io, M)
  io = pretty(io)
  print(io, "Map defined by a Julia function")
  if isdefined(M, :g)
    print(io, " with inverse")
  end
  println(io)
  println(io, Indent(),"from ", Lowercase(), domain(M))
  print(io, "to ", Lowercase(), codomain(M), Dedent())
end

function MapFromFunc(D, C, f)
  return MapFromFunc{typeof(D), typeof(C)}(D, C, f)
end

function MapFromFunc(D, C, f, g)
  return MapFromFunc{typeof(D), typeof(C)}(D, C, f, g)
end

function Base.inv(M::MapFromFunc)
  if isdefined(M, :g)
     return MapFromFunc(codomain(M), domain(M), M.g, M.f)
  else
     return MapFromFunc(codomain(M), domain(M), x->preimage(M, x))
  end
end

###########################################################
# InverseMap
###########################################################

mutable struct InverseMap{D, C} <: Map{D, C, MapWithHeader, InverseMap}
  header::MapHeader{D, C}
  origin::Map{C, D}

  function InverseMap{D, C}(f::Map{C, D}) where {D, C}
    z = new{D, C}()
    z.header = MapHeader(codomain(f), domain(f), inverse_fn(f), image_fn(f))
    z.origin = f
    return z
  end
end

function InverseMap(f::Map{C, D}) where {D, C}
  return InverseMap{D, C}(f)
end

function show(io::IO, M::InverseMap)
  @show_name(io, M)
  println(io, "inverse of")
  print(io, " ", M.origin)
end
