###############################################################################
#
#   PuiseuxSeries.jl : Puiseux series over rings and fields
#
###############################################################################

###############################################################################
#
#   PuiseuxSeriesRing constructor
#
###############################################################################

@doc raw"""
    puiseux_series_ring(R::Ring, prec::Int, s::VarName; cached::Bool=true)
    puiseux_series_field(R::Field, prec::Int, s::VarName; cached::Bool=true)

Return a tuple $(S, x)$ consisting of the parent object `S` of a Puiseux series
ring over the given base ring and a generator `x` for the Puiseux series ring.
The maximum precision of the series in the ring is set to `prec`. This is taken as a
maximum relative precision of the underlying Laurent series that are used to implement
the Puiseux series in the ring. The supplied string `s` specifies the way the
generator of the Puiseux series ring will be printed. By default, the parent
object `S` will be cached so that supplying the same base ring, string and
precision in future will return the same parent object and generator. If
caching of the parent object is not required, `cached` can be set to `false`.
"""
puiseux_series_ring(R::Ring, prec::Int, s::VarName; cached::Bool=true) =
   Generic.PuiseuxSeriesRing(R, prec, Symbol(s); cached)

puiseux_series_field(R::Field, prec::Int, s::VarName; cached::Bool=true) =
   Generic.PuiseuxSeriesField(R, prec, Symbol(s); cached)
