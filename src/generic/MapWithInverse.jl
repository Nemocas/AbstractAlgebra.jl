################################################################################
#
#  MapWithInverse.jl : Map with section, retraction, two-sided inverse, etc.
#
################################################################################

export map_with_preimage_from_func, map_with_section_from_func, 
       map_with_retraction_from_func, image_map, preimage_map, section_map,
       retraction_map, map_with_retraction, map_with_section, inverse_fn

################################################################################
#
#  MapWithSection
#
################################################################################

map_with_section(f::Map{D, C}, g::Map{C, D}) where {D, C} = MapWithSection(f, g)

domain(f::MapWithSection{D, C}) where {D, C} = get_field(f.map, :domain)::D
codomain(f::MapWithSection{D, C}) where {D, C} = get_field(f.map, :codomain)::C
image_fn(f::MapWithSection) = image_fn(f.map)
inverse_fn(f::MapWithSection) = image_fn(f.section)
image_map(f::MapWithSection) = f.map
preimage_map(f::MapWithSection) = f.section # for convenience only
section_map(f::MapWithSection) = f.section

(f::MapWithSection{D, C})(a) where {D, C} = (f.map)(a)::elem_type(C)

# These two functions are provided for convenience only. Strictly speaking
# preimage is not the correct name for this type of construction.
function map_with_preimage_from_func(domain, codomain, image_fn::Function, inverse_fn::Function)
   return MapWithSection(FunctionalMap(domain, codomain, image_fn),
                          FunctionalMap(codomain, domain, inverse_fn))
end

function map_with_preimage_from_func(domain, codomain, image_fn::Function)
   return MapWithSection(FunctionalMap(domain, codomain, image_fn))
end

function map_with_section_from_func(domain, codomain, image_fn::Function, inverse_fn::Function)
   return MapWithSection(FunctionalMap(domain, codomain, image_fn),
                          FunctionalMap(codomain, domain, inverse_fn))
end

function map_with_section_from_func(domain, codomain, image_fn::Function)
   return MapWithSection(FunctionalMap(domain, codomain, image_fn))
end

function show(io::IO, M::MapWithSection)
   println(io, "Map with section with the following data")
   println(io, "")
   println(io, "Domain:")
   println(io, "=======")
   println(io, domain(M))
   println(io, "")
   println(io, "Codomain:")
   println(io, "========")
   println(io, codomain(M))
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

function inv(f::MapWithSection)
   return MapWithRetraction(f.section, f.map)
end
 
################################################################################
#
#  MapWithRetraction
#
################################################################################

map_with_retraction(f::Map{D, C}, g::Map{C, D}) where {D, C} = MapWithRetraction(f, g)

domain(f::MapWithRetraction{D, C}) where {D, C} = get_field(f.map, :domain)::D
codomain(f::MapWithRetraction{D, C}) where {D, C} = get_field(f.map, :codomain)::C
image_fn(f::MapWithRetraction) = image_fn(f.map)
inverse_fn(f::MapWithRetraction) = image_fn(f.retraction)
image_map(f::MapWithRetraction) = f.map
retraction_map(f::MapWithRetraction) = f.retraction

retraction_map(f::MapCache) = retraction_map(f.map)

(f::MapWithRetraction{D, C})(a) where {D, C} = (f.map)(a)::elem_type(C)

function map_with_retraction_from_func(domain, codomain, image_fn::Function, inverse_fn::Function)
   return MapWithRetraction(FunctionalMap(domain, codomain, image_fn),
                          FunctionalMap(codomain, domain, inverse_fn))
end

function map_with_retraction_from_func(domain, codomain, image_fn::Function)
   return MapWithRetraction(FunctionalMap(domain, codomain, image_fn))
end

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

function inv(f::MapWithRetraction)
   return MapWithSection(f.retraction, f.map)
end

