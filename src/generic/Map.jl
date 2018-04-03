################################################################################
#
#  Map.jl : Maps
#
################################################################################

export MapFromFunc, compose, domain, codomain, image

function check_composable(a::AbstractAlgebra.Map{U, C}, b::AbstractAlgebra.Map{D, U}) where {C, U, D}
   domain(a) == codomain(b) || throw("Incompatible maps")
end

################################################################################
#
#  IdentityMap
#
################################################################################

domain(f::IdentityMap) = f.domain
codomain(f::IdentityMap) = f.codomain
image(f::IdentityMap) = x -> x

(f::IdentityMap)(a) = a

################################################################################
#
#  FunctionalMap
#
################################################################################

domain(f::FunctionalMap) = f.domain
codomain(f::FunctionalMap) = f.codomain
image(f::FunctionalMap) = f.image

(f::FunctionalMap)(a) = image(f)(a)

function MapFromFunc(image::Function, domain::D, codomain::C) where {C, D}
   return FunctionalMap(domain, codomain, image)
end

################################################################################
#
#  CompositeMap
#
################################################################################

domain(f::CompositeMap) = domain(f.map1)
codomain(f::CompositeMap) = codomain(f.map2)

# This is a device to prevent Julia trying to compute
# the types of a very long composition of closures
struct UntypedFunction <: Function
   fn::Function
end

function compose(f::UntypedFunction, g::UntypedFunction)
   return x -> f.fn(g.fn(x))
end

(f::UntypedFunction)(a) = f.fn(a)

function image(f::CompositeMap)
   if isdefined(f, :fn_cache)
      return f.fn_cache
   else
      f1 = image(f.map1)
      f2 = image(f.map2)
      fn = compose(UntypedFunction(f2), UntypedFunction(f1))
      f.fn_cache = fn
      return fn
   end
end

(f::CompositeMap)(a) = image(f)(a)

function compose(f::AbstractAlgebra.Map{U, C}, g::AbstractAlgebra.Map{D, U}) where {D, U, C}
   check_composable(f, g)
   return CompositeMap(f, g)
end

function compose(f::FunctionalMap{U, C}, g::FunctionalMap{D, U}) where {D, U,
C}
   check_composable(f, g)
   return CompositeMap(f, g)
end

function compose(f::AbstractAlgebra.Map{D, C}, g::AbstractAlgebra.IdentityMap{D}) where {D, C}
   check_composable(f, g)
   return f
end

function compose(f::AbstractAlgebra.IdentityMap{C}, g::AbstractAlgebra.Map{D, C}) where {D, C}
   check_composable(f, g)
   return g
end


