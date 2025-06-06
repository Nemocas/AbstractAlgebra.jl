################################################################################
#
#  MapWithInverse.jl : Generic Map with section, retraction,
#                      two-sided inverse, etc.
#
################################################################################

################################################################################
#
#  MapWithSection
#
################################################################################

domain(f::MapWithSection{D, C}) where {D, C} = domain(f.map)::D
codomain(f::MapWithSection{D, C}) where {D, C} = codomain(f.map)::C
image_fn(f::MapWithSection) = image_fn(f.map)
inverse_fn(f::MapWithSection) = image_fn(f.section)
image_map(f::MapWithSection) = f.map
preimage_map(f::MapWithSection) = f.section # for convenience only
section_map(f::MapWithSection) = f.section

(f::MapWithSection{D, C})(a) where {D, C} = image(f, a)::elem_type(C)

function preimage(f::MapWithSection{D, C}, a) where {D, C}
  return inverse_fn(f)(a)::elem_type(D)
end

function image(f::MapWithSection{D, C}, a) where {D, C}
  return image_fn(f)(a)::elem_type(C)
end

function Base.show(io::IO, M::MapWithSection)
   if is_terse(io)
      print(io, "Map with section")
   else
      io = pretty(io)
      io = terse(io)
      print(io, "Map: ")
      print(io, Lowercase(), domain(M), " -> ")
      print(io, Lowercase(), codomain(M))
   end
end

function compose(f::MapWithSection{U, C}, g::MapWithSection{D, U}) where {D, U, C}
   check_composable(f, g)
   m = compose(f.map, g.map)
   if isdefined(g, :section) && isdefined(f, :section)
      p = compose(g.section, f.section)
      return MapWithSection(m, p)
   else
      return MapWithSection(m)
   end
end

function Base.inv(f::MapWithSection)
   return MapWithRetraction(f.section, f.map)
end

################################################################################
#
#  MapWithRetraction
#
################################################################################

domain(f::MapWithRetraction{D, C}) where {D, C} = domain(f.map)::D
codomain(f::MapWithRetraction{D, C}) where {D, C} = codomain(f.map)::C
image_fn(f::MapWithRetraction) = image_fn(f.map)
inverse_fn(f::MapWithRetraction) = image_fn(f.retraction)
image_map(f::MapWithRetraction) = f.map
retraction_map(f::MapWithRetraction) = f.retraction

retraction_map(f::MapCache) = retraction_map(f.map)

(f::MapWithRetraction{D, C})(a) where {D, C} = image(f, a)::elem_type(C)

function preimage(f::MapWithRetraction{D, C}, a) where {D, C}
  return inverse_fn(f)(a)::elem_type(D)
end

function image(f::MapWithRetraction{D, C}, a) where {D, C}
  return image_fn(f)(a)::elem_type(C)
end

function Base.show(io::IO, M::MapWithRetraction)
   if is_terse(io)
      print(io, "Map with retraction")
   else
      io = pretty(io)
      io = terse(io)
      print(io, "Map: ")
      print(io, Lowercase(), domain(M), " -> ")
      print(io, Lowercase(), codomain(M))
   end
end

function compose(f::MapWithRetraction{U, C}, g::MapWithRetraction{D, U}) where {D, U, C}
   check_composable(f, g)
   m = compose(f.map, g.map)
   if isdefined(g, :retraction) && isdefined(f, :retraction)
      p = compose(g.retraction, f.retraction)
      return MapWithRetraction(m, p)
   else
      return MapWithRetraction(m)
   end
end

function Base.inv(f::MapWithRetraction)
   return MapWithSection(f.retraction, f.map)
end
