################################################################################
#
#  MapWith.jl : Map with preimage, left inverse, right inverse, etc.
#
################################################################################

export MapWithPreimageFromFunc, preimage_map, image_map

################################################################################
#
#  MapWithPreimage
#
################################################################################

domain(f::MapWithPreimage) = domain(f.image_map)
codomain(f::MapWithPreimage) = codomain(f.image_map)
image_fn(f::MapWithPreimage) = image_fn(f.image_map)
preimage_fn(f::MapWithPreimage) = image_fn(f.preimage_map)
image_map(f::MapWithPreimage) = f.image_map
preimage_map(f::MapWithPreimage) = f.preimage_map

(f::MapWithPreimage)(a) = image_fn(f.image_map)(a)

function MapWithPreimageFromFunc(image_fn::Function, preimage_fn::Function, domain, codomain)
   return MapWithPreimage(FunctionalMap(domain, codomain, image_fn),
                          FunctionalMap(codomain, domain, preimage_fn))
end

function MapWithPreimageFromFunc(image_fn::Function, domain, codomain)
   return MapWithPreimage(FunctionalMap(domain, codomain, image_fn))
end

function show(io::IO, M::MapWithPreimage)
   println(io, "Map with preimage with the following data")
   println(io, "")
   println(io, "Domain:")
   println(io, "=======")
   println(io, domain(M))
   println(io, "")
   println(io, "Coomain:")
   println(io, "========")
   println(io, codomain(M))
end

function compose(f::MapWithPreimage{U, C}, g::MapWithPreimage{D, U}) where {D, U, C}
   check_composable(f, g)
   m = compose(f.image_map, g.image_map)
   if isdefined(g, :preimage_map) && isdefined(f, :preimage_map)
      p = compose(g.preimage_map, f.preimage_map)
      return MapWithPreimage(m, p)
   else
      return MapWithPreimage(m)
   end
end

