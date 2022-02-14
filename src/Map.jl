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

domain(f::Map) = f.domain # fallback, FIXME: to be removed eventually
codomain(f::Map) = f.codomain # fallback, FIXME: to be removed eventually
image_fn(f::Map) = f.image_fn # fallback, FIXME: to be removed eventually


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

