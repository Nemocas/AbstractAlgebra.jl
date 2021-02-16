##############################################################################
#
#  Parent object constructors
#
##############################################################################

function Nemo.PolynomialRing(R::Nemo.Ring, n::Int, s::String="x";
    cached::Bool = false, ordering::Symbol = :lex)
    return Nemo.PolynomialRing(R, ["$s$i" for i=1:n], cached = cached,
                           ordering = ordering)
end