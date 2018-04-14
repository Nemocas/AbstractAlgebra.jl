################################################################################
#
#  Map.jl : Maps
#
################################################################################

export map_from_func, compose, domain, codomain

################################################################################
#
#  All maps
#
################################################################################

function check_composable(a::AbstractAlgebra.Map{U, C}, b::AbstractAlgebra.Map{D, U}) where {C, U, D}
   domain(a) != codomain(b) && error("Incompatible maps")
end

###############################################################################
#
#   CompositeMap
#
###############################################################################

domain(f::CompositeMap) = f.domain
codomain(f::CompositeMap) = f.codomain

function (f::CompositeMap)(a)
   return f.map2(f.map1(a))
end

function compose(f::AbstractAlgebra.Map{U, C}, g::AbstractAlgebra.Map{D, U}) where {D, U, C}
   check_composable(f, g)
   return CompositeMap(f, g)
end

function show_short(io::IO, M::CompositeMap)
   show_short(io, M.map2)
   println(io, "then")
   show(io, M.map1)
end

function show(io::IO, M::CompositeMap)
   println(io, "Composite map consisting of the following")
   println(io, "")
   show_short(io, M.map2)
   println(io, "then")
   show_short(io, M.map1)
end

################################################################################
#
#  IdentityMap
#
################################################################################

domain(f::IdentityMap) = f.domain
codomain(f::IdentityMap) = f.codomain
image_fn(f::IdentityMap) = x -> x

(f::IdentityMap)(a) = a

function show(io::IO, M::IdentityMap)
   println(io, "Identity map with")
   println(io, "")
   println(io, "Domain:")
   println(io, "=======")
   println(io, domain(M))
end

function compose(f::AbstractAlgebra.Map{C, C, <:AbstractAlgebra.IdentityMap}, g::AbstractAlgebra.Map{D, C}) where {D, C}
   check_composable(f, g)
   return g
end

function compose(f::AbstractAlgebra.Map{D, C}, g::AbstractAlgebra.Map{D, D, <:AbstractAlgebra.IdentityMap}) where {D, C}
   check_composable(f, g)
   return f
end

function compose(f::AbstractAlgebra.Map{D, D, <:AbstractAlgebra.IdentityMap}, g::AbstractAlgebra.Map{D, D, <:AbstractAlgebra.IdentityMap}) where D
   check_composable(f, g)
   return g
end

################################################################################
#
#  FunctionalMap
#
################################################################################

domain(f::FunctionalMap) = f.domain
codomain(f::FunctionalMap) = f.codomain
image_fn(f::FunctionalMap) = f.image_fn

image_fn(f::MapCache) = image_fn(f.map)

function (f::FunctionalMap)(a)
   parent(a) != domain(f) && throw(DomainError())
   return image_fn(f)(a)
end

function map_from_func(domain, codomain, image_fn::Function)
   return FunctionalMap(domain, codomain, image_fn)
end

function show(io::IO, M::FunctionalMap)
   println(io, "Map with the following data")
   println(io, "")
   println(io, "Domain:")
   println(io, "=======")
   println(io, domain(M))
   println(io, "")
   println(io, "Codomain:")
   println(io, "========")
   println(io, codomain(M))
end

################################################################################
#
#  FunctionalCompositeMap
#
################################################################################

domain(f::FunctionalCompositeMap) = f.domain
codomain(f::FunctionalCompositeMap) = f.codomain

# This is a device to prevent Julia trying to compute
# the types of a very long composition of closures
struct UntypedFunction <: Function
   fn::Function
end

function compose(f::UntypedFunction, g::UntypedFunction, d)
   return function(x)
      parent(x) != d && error("Element not in domain of map")
      return f.fn(g.fn(x))
   end
end

(f::UntypedFunction)(a) = f.fn(a)

function image_fn(f::FunctionalCompositeMap)
   if isdefined(f, :fn_cache)
      return f.fn_cache
   else
      f1 = image_fn(f.map1)
      f2 = image_fn(f.map2)
      d = domain(f)
      fn = compose(UntypedFunction(f2), UntypedFunction(f1), d)
      f.fn_cache = fn
      return fn
   end
end

function (f::FunctionalCompositeMap)(a)
   return image_fn(f)(a)
end

function compose(f::Map{U, C, <:AbstractAlgebra.FunctionalMap}, g::Map{D, U, <:AbstractAlgebra.FunctionalMap}) where {D, U, C}
   check_composable(f, g)
   return FunctionalCompositeMap(f, g)
end

function show_short(io::IO, M::AbstractAlgebra.Map)
   println(domain(M), " -> ", codomain(M))
end

function show_short(io::IO, M::FunctionalCompositeMap)
   show_short(io, M.map2)
   println(io, "then")
   show(io, M.map1)
end

function show(io::IO, M::FunctionalCompositeMap)
   println(io, "Composite map consisting of the following")
   println(io, "")
   show_short(io, M.map2)
   println(io, "then")
   show_short(io, M.map1)
end

