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

get_field(M, f) = getfield(M, f) # fall back to Julia builtin
set_field!(M, f, g) = setfield!(M, f, g) # fall back to Julia builtin

domain(f::Map) = get_field(f, :domain)
codomain(f::Map) = get_field(f, :codomain)

function check_composable(a::Map, b::Map)
   codomain(a) != domain(b) && error("Incompatible maps")
end

###############################################################################
#
#   CompositeMap
#
###############################################################################

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

identity_map(R::D) where D <: AbstractAlgebra.Set = Generic.IdentityMap{D}(R)

################################################################################
#
#  FunctionalMap
#
################################################################################

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

