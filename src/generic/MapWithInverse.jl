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

(f::MapWithSection{D, C})(a) where {D, C} = (f.map)(a)::elem_type(C)

function show(io::IO, M::MapWithSection)
   println(io, "Map with section with the following data")
   println(io, "")
   println(io, "Domain:")
   println(io, "=======")
   println(io, domain(M))
   println(io, "")
   println(io, "Codomain:")
   println(io, "========")
   print(io, codomain(M))
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

(f::MapWithRetraction{D, C})(a) where {D, C} = (f.map)(a)::elem_type(C)

function show(io::IO, M::MapWithRetraction)
   println(io, "Map with retraction with the following data")
   println(io, "")
   println(io, "Domain:")
   println(io, "=======")
   println(io, domain(M))
   println(io, "")
   println(io, "Codomain:")
   println(io, "========")
   println(io, codomain(M))
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
