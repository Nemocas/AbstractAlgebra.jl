################################################################################
#
#  MapWithInverse.jl : Map with section, retraction, two-sided inverse, etc.
#
################################################################################

################################################################################
#
#  MapWithSection
#
################################################################################

@doc raw"""
    map_with_section(f::Map{D, C}, s::Map{C, D}) where {D, C}

Return the map `f` together with a known section `s` of it, i.e. a map with
$f(s(y)) = y$ for all $y$ in the codomain of `f`.
"""
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

@doc raw"""
    map_with_section_from_func(f::Function, s::Function, R, S)

Return the map from `R` to `S` given by the Julia function `f`, together with
the section given by the Julia function `s`. See [`map_with_section`](@ref).
"""
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

@doc raw"""
    map_with_retraction(f::Map{D, C}, r::Map{C, D}) where {D, C}

Return the map `f` together with a known retraction `r` of it, i.e. a map with
$r(f(x)) = x$ for all $x$ in the domain of `f`.
"""
map_with_retraction(f::Map{D, C}, g::Map{C, D}) where {D, C} = Generic.MapWithRetraction(f, g)

@doc raw"""
    map_with_retraction_from_func(f::Function, r::Function, R, S)

Return the map from `R` to `S` given by the Julia function `f`, together with
the retraction given by the Julia function `r`. See
[`map_with_retraction`](@ref).

# Examples

```jldoctest
julia> f = map_with_retraction_from_func(x -> x + 1, x -> x - 1, ZZ, ZZ)
Map with retraction
  from integers
  to integers

julia> f(ZZ(1))
2
```
"""
function map_with_retraction_from_func(image_fn::Function, inverse_fn::Function, domain, codomain)
   return Generic.MapWithRetraction(Generic.FunctionalMap(domain, codomain, image_fn),
                          Generic.FunctionalMap(codomain, domain, inverse_fn))
end

function map_with_retraction_from_func(image_fn::Function, domain, codomain)
   return Generic.MapWithRetraction(Generic.FunctionalMap(domain, codomain, image_fn))
end

