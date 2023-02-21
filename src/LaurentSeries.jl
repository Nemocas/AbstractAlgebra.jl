###############################################################################
#
#   LaurentSeries.jl : Laurent series over rings and fields,
#                      capped relative precision
#
###############################################################################

export laurent_series_field, laurent_series_ring, deflate, downscale, exp_gcd,
       inflate, upscale

###############################################################################
#
#   power_series_ring constructor
#
###############################################################################

@doc Markdown.doc"""
    laurent_series_ring(R::Ring, prec::Int, s::Union{Char, AbstractString, Symbol}; cached=true)
    laurent_series_field(R::Field, prec::Int, s::Union{Char, AbstractString; cached=true)

Return a tuple $(S, x)$ consisting of the parent object `S` of a Laurent series
ring over the given base ring and a generator `x` for the Laurent series ring.
The maximum precision of the series in the ring is set to `prec`. This is taken as a
maximum relative precision. The supplied string `s` specifies the way the
generator of the Laurent series ring will be printed. By default, the parent
object `S` will be cached so that supplying the same base ring, string and
precision in future will return the same parent object and generator. If
caching of the parent object is not required, `cached` can be set to `false`.
"""
function laurent_series_ring(R::Ring, prec::Int, s::AbstractString; cached=true)
   return laurent_series_ring(R, prec, Symbol(s); cached=cached)
end

function laurent_series_ring(R::Ring, prec::Int, s::Symbol; cached=true)
   return Generic.laurent_series_ring(R, prec, s; cached=cached)
end

function laurent_series_ring(R::Ring, prec::Int, s::Char; cached=true)
   return laurent_series_ring(R, prec, Symbol(s); cached=cached)
end

function laurent_series_field(R::Field, prec::Int, s::AbstractString; cached=true)
   return laurent_series_field(R, prec, Symbol(s); cached=cached)
end

function laurent_series_field(R::Field, prec::Int, s::Symbol; cached=true)
   return Generic.laurent_series_field(R, prec, s; cached=cached)
end

function laurent_series_field(R::Field, prec::Int, s::Char; cached=true)
   return laurent_series_field(R, prec, Symbol(s); cached=cached)
end

