###############################################################################
#
#   Random generation
#
###############################################################################

rand(x::Union{Integers,Rationals,Floats}, v) = rand(Random.GLOBAL_RNG, x, v)

rand(x::Union{AbstractAlgebra.MatSpace,AbstractAlgebra.MatAlgebra,AbstractAlgebra.FPModule,
              AbstractAlgebra.ResField,AbstractAlgebra.FracField,AbstractAlgebra.ResRing},
     v...) = rand(Random.GLOBAL_RNG, x, v...)

rand(x::Union{SeriesRing,AbstractAlgebra.NCPolyRing,AbstractAlgebra.PolyRing,
              Generic.LaurentSeriesRing,Generic.LaurentSeriesField},
     v1, v...) = rand(Random.GLOBAL_RNG, x, v1, v...)

rand(x::Union{Generic.PuiseuxSeriesRing,Generic.PuiseuxSeriesField,AbstractAlgebra.MPolyRing}, v1, v2, v...) =
   rand(Random.GLOBAL_RNG, x, v1, v2, v...)
