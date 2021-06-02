###############################################################################
#
#   PuiseuxSeries.jl : Puiseux series over rings and fields
#
###############################################################################

export PuiseuxSeriesField, PuiseuxSeriesRing, laurent_ring, rescale!

###############################################################################
#
#   PuiseuxSeriesRing constructor
#
###############################################################################

@doc Markdown.doc"""
    PuiseuxSeriesRing(R::Ring, prec::Int, s::Union{Symbol, Char, AbstractString}; cached=true)
    PuiseuxSeriesField(R::Field, prec::Int, s::Union{Symbol, Char, AbstractString}; cached = true)

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
function PuiseuxSeriesRing(R::Ring, prec::Int, s::AbstractString; cached=true)
   return Generic.PuiseuxSeriesRing(R, prec, Symbol(s); cached=cached)
end

function PuiseuxSeriesRing(R::Ring, prec::Int, s::Symbol; cached=true)
   return Generic.PuiseuxSeriesRing(R, prec, s; cached=cached)
end

function PuiseuxSeriesRing(R::Ring, prec::Int, s::Char; cached=true)
   return Generic.PuiseuxSeriesRing(R, prec, Symbol(s); cached=cached)
end

function PuiseuxSeriesField(R::Field, prec::Int, s::AbstractString; cached = true)
   return Generic.PuiseuxSeriesField(R, prec, Symbol(s); cached=cached)
end

function PuiseuxSeriesField(R::Field, prec::Int, s::Symbol; cached = true)
   return Generic.PuiseuxSeriesField(R, prec, s; cached=cached)
end

function PuiseuxSeriesField(R::Field, prec::Int, s::Char; cached = true)
   return Generic.PuiseuxSeriesField(R, prec, Symbol(s); cached=cached)
end
