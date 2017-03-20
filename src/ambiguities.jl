# Fix ambiguities on julia 0.4

*(a::ResElem{fmpz}, b::fmpz) = parent(a)(data(a) * b)

*(a::fmpz, b::ResElem{fmpz}) = b*a

+(a::ResElem{fmpz}, b::fmpz) = parent(a)(data(a) + b)

+(a::fmpz, b::ResElem{fmpz}) = b + a

-(a::ResElem{fmpz}, b::fmpz) = parent(a)(data(a) - b)

-(a::fmpz, b::ResElem{fmpz}) = parent(b)(a - data(b))

function ==(a::ResElem{fmpz}, b::fmpz)
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

==(a::fmpz, b::ResElem{fmpz}) = b == a

#

*(::fmpz, ::PolyElem{fmpz}) = nothing

*(::PolyElem{fmpz}, ::fmpz) = nothing

+(::fmpz, ::PolyElem{fmpz}) = nothing

+(::PolyElem{fmpz}, ::fmpz) = nothing

-(::fmpz, ::PolyElem{fmpz}) = nothing

-(::PolyElem{fmpz}, ::fmpz) = nothing

==(::fmpz, ::PolyElem{fmpz}) = nothing

==(::PolyElem{fmpz}, ::fmpz) = nothing

divexact(::PolyElem{fmpz}, ::fmpz) = nothing

evaluate(::PolyElem{fmpz}, ::fmpz) = nothing

#

*(::fmpz, ::SeriesElem{fmpz}) = nothing

*(::SeriesElem{fmpz}, ::fmpz) = nothing

+(::fmpz, ::SeriesElem{fmpz}) = nothing

+(::SeriesElem{fmpz}, ::fmpz) = nothing

-(::fmpz, ::SeriesElem{fmpz}) = nothing

-(::SeriesElem{fmpz}, ::fmpz) = nothing

==(::fmpz, ::SeriesElem{fmpz}) = nothing

==(::SeriesElem{fmpz}, ::fmpz) = nothing

*(::fmpz, ::RelSeriesElem{fmpz}) = nothing

*(::RelSeriesElem{fmpz}, ::fmpz) = nothing

+(::fmpz, ::RelSeriesElem{fmpz}) = nothing

+(::RelSeriesElem{fmpz}, ::fmpz) = nothing

-(::fmpz, ::RelSeriesElem{fmpz}) = nothing

-(::RelSeriesElem{fmpz}, ::fmpz) = nothing

==(::fmpz, ::RelSeriesElem{fmpz}) = nothing

==(::RelSeriesElem{fmpz}, ::fmpz) = nothing

*(::fmpz, ::AbsSeriesElem{fmpz}) = nothing

*(::AbsSeriesElem{fmpz}, ::fmpz) = nothing

+(::fmpz, ::AbsSeriesElem{fmpz}) = nothing

+(::AbsSeriesElem{fmpz}, ::fmpz) = nothing

-(::fmpz, ::AbsSeriesElem{fmpz}) = nothing

-(::AbsSeriesElem{fmpz}, ::fmpz) = nothing

==(::fmpz, ::AbsSeriesElem{fmpz}) = nothing

==(::AbsSeriesElem{fmpz}, ::fmpz) = nothing

*(::fmpz, ::MatElem{fmpz}) = nothing

*(::MatElem{fmpz}, ::fmpz) = nothing

+(::fmpz, ::MatElem{fmpz}) = nothing

+(::MatElem{fmpz}, ::fmpz) = nothing

-(::fmpz, ::MatElem{fmpz}) = nothing

-(::MatElem{fmpz}, ::fmpz) = nothing

==(::MatElem{fmpz}, ::fmpz) = nothing

divexact(::MatElem{fmpz}, ::fmpz) = nothing

#

setindex_t!(a::nmod_mat, b::GenRes{fmpz}, i::Int, j::Int) = setindex_!(a, data(b), i, j)

*(::FracElem{fmpz}, ::fmpz) = nothing

*(::fmpz, ::FracElem{fmpz}) = nothing

+(::FracElem{fmpz}, ::fmpz) = nothing

+(::fmpz, ::FracElem{fmpz}) = nothing

-(::FracElem{fmpz}, ::fmpz) = nothing

-(::fmpz, ::FracElem{fmpz}) = nothing

==(::FracElem{fmpz}, ::fmpz) = nothing

==(::fmpz, ::FracElem{fmpz}) = nothing

divexact(::FracElem{fmpz}, ::fmpz) = nothing

divexact(::fmpz, ::FracElem{fmpz}) = nothing

