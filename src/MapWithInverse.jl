################################################################################
#
#  MapWithInverse.jl : Map with section, retraction, two-sided inverse, etc.
#
################################################################################

export image_map, inverse_fn, map_with_preimage_from_func, map_with_retraction,
       map_with_retraction_from_func, map_with_section,
       map_with_section_from_func, preimage_map, retraction_map, section_map
       
################################################################################
#
#  MapWithSection
#
################################################################################

map_with_section(f::Map{D, C}, g::Map{C, D}) where {D, C} = Generic.MapWithSection(f, g)

# These two functions are provided for convenience only. Strictly speaking
# preimage is not the correct name for this type of construction.
function map_with_preimage_from_func(image_fn::Function, inverse_fn::Function, domain, codomain)
   return Generic.MapWithSection(Generic.FunctionalMap(domain, codomain, image_fn),
                          Generic.FunctionalMap(codomain, domain, inverse_fn))
end

function map_with_preimage_from_func(image_fn::Function, domain, codomain)
   return Generic.MapWithSection(Generic.FunctionalMap(domain, codomain, image_fn))
end

function map_with_section_from_func(image_fn::Function, inverse_fn::Function, domain, codomain)
   return Generic.MapWithSection(Generic.FunctionalMap(domain, codomain, image_fn),
                          Generic.FunctionalMap(codomain, domain, inverse_fn))
end

function map_with_section_from_func(image_fn::Function, domain, codomain)
   return Generic.MapWithSection(Generic.FunctionalMap(domain, codomain, image_fn))
end
 
################################################################################
#
#  MapWithRetraction
#
################################################################################

map_with_retraction(f::Map{D, C}, g::Map{C, D}) where {D, C} = Generic.MapWithRetraction(f, g)

function map_with_retraction_from_func(image_fn::Function, inverse_fn::Function, domain, codomain)
   return Generic.MapWithRetraction(Generic.FunctionalMap(domain, codomain, image_fn),
                          Generic.FunctionalMap(codomain, domain, inverse_fn))
end

function map_with_retraction_from_func(image_fn::Function, domain, codomain)
   return Generic.MapWithRetraction(Generic.FunctionalMap(domain, codomain, image_fn))
end

