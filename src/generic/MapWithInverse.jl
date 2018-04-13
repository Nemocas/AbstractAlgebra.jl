################################################################################
#
#  MapWithInverse.jl : Map with section, retraction, two-sided inverse, etc.
#
################################################################################

export map_with_preimage_from_func, map_with_section_from_func, 
       map_with_retraction_from_func, image_map, preimage_map, section_map,
       retraction_map

################################################################################
#
#  MapWithSection
#
################################################################################

domain(f::MapWithSection) = domain(f.map)
codomain(f::MapWithSection) = codomain(f.map)
image_fn(f::MapWithSection) = image_fn(f.map)
inverse_fn(f::MapWithSection) = image_fn(f.section)
image_map(f::MapWithSection) = f.map
preimage_map(f::MapWithSection) = f.section # for convenience only
section_map(f::MapWithSection) = f.section

inverse_fn(f::MapCache) = inverse_fn(f.map)
inverse_map(f::MapCache) = inverse_map(f.map)
preimage_map(f::MapCache) = preimage_map(f.map)
section_map(f::MapCache) = section_map(f.map)

(f::MapWithSection)(a) = (f.map)(a)

# These two functions are provided for convenience only. Strictly speaking
# preimage is not the correct name for this type of construction.
function map_with_preimage_from_func(image_fn::Function, inverse_fn::Function, domain, codomain)
   return MapWithSection(FunctionalMap(image_fn, domain, codomain),
                          FunctionalMap(inverse_fn, codomain, domain))
end

function map_with_preimage_from_func(image_fn::Function, domain, codomain)
   return MapWithSection(FunctionalMap(image_fn, domain, codomain))
end

function map_with_section_from_func(image_fn::Function, inverse_fn::Function, domain, codomain)
   return MapWithSection(FunctionalMap(image_fn, domain, codomain),
                          FunctionalMap(inverse_fn, codomain, domain))
end

function map_with_section_from_func(image_fn::Function, domain, codomain)
   return MapWithSection(FunctionalMap(image_fn, domain, codomain))
end

function show(io::IO, M::MapWithSection)
   println(io, "Map with section with the following data")
   println(io, "")
   println(io, "Domain:")
   println(io, "=======")
   println(io, domain(M))
   println(io, "")
   println(io, "Coomain:")
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

################################################################################
#
#  MapWithRetraction
#
################################################################################

domain(f::MapWithRetraction) = domain(f.map)
codomain(f::MapWithRetraction) = codomain(f.map)
image_fn(f::MapWithRetraction) = image_fn(f.map)
inverse_fn(f::MapWithRetraction) = image_fn(f.retraction)
image_map(f::MapWithRetraction) = f.map
retraction_map(f::MapWithSection) = f.retraction

retraction_map(f::MapCache) = retraction_map(f.map)

(f::MapWithRetraction)(a) = (f.map)(a)

function map_with_retraction_from_func(image_fn::Function, inverse_fn::Function, domain, codomain)
   return MapWithRetraction(FunctionalMap(image_fn, domain, codomain),
                          FunctionalMap(inverse_fn, codomain, domain))
end

function map_with_retraction_from_func(image_fn::Function, domain, codomain)
   return MapWithRetraction(FunctionalMap(image_fn, domain, codomain))
end

function show(io::IO, M::MapWithRetraction)
   println(io, "Map with retraction with the following data")
   println(io, "")
   println(io, "Domain:")
   println(io, "=======")
   println(io, domain(M))
   println(io, "")
   println(io, "Coomain:")
   println(io, "========")
   println(io, codomain(M))
end

