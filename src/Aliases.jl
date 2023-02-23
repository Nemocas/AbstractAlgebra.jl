# for backwards compatibility
@alias iscompatible is_compatible
@alias isconstant is_constant
@alias isdegree is_degree
@alias isdivisible_by is_divisible_by
@alias isdomain_type is_domain_type
@alias isexact_type is_exact_type
@alias isgen is_gen
@alias ishessenberg is_hessenberg
@alias ishnf is_hnf
@alias ishomogeneous is_homogeneous
@alias isinvertible is_invertible
@alias isinvertible_with_inverse is_invertible_with_inverse
@alias isirreducible is_irreducible
@alias isisomorphic is_isomorphic
@alias ismonic is_monic
@alias ismonomial is_monomial
@alias ismonomial_recursive is_monomial_recursive
@alias isnegative is_negative
@alias ispopov is_popov
@alias ispower is_power
@alias ispower_with_root is_power_with_root
@alias isprobable_prime is_probable_prime
@alias isreverse is_reverse
@alias isrimhook is_rimhook
@alias isrref is_rref
@alias issnf is_snf
@alias issquare is_square
@alias issquare_with_sqrt is_square_with_sqrt
@alias issquarefree is_squarefree
@alias issubmodule is_submodule
@alias isterm is_term
@alias isterm_recursive is_term_recursive
@alias isunit is_unit
@alias isunivariate is_univariate
@alias isweak_popov is_weak_popov
@alias iszero_column is_zero_column
@alias iszero_row is_zero_row
@alias Localization localization

#
# Some deprecation aliases follow. They use @alias instead of @deprecate
# because there are still packages adding methods for them, which would break
# if we used @deprecate. Once those packages are updated, we can switch to
# using @deprecate
#

# Deprecated in 0.27.*

@alias MatrixSpace matrix_space
#@deprecate MatrixSpace(R::Ring, rows::Int, cols::Int; cache::Bool=true) matrix_space(R, rows, cols; cache)

@alias PolynomialRing polynomial_ring
@alias PowerSeriesRing power_series_ring
@alias LaurentSeriesRing laurent_series_ring
@alias LaurentSeriesField laurent_series_field
@alias FractionField fraction_field
@alias NumberField number_field

@alias PolyElem PolyRingElem
@alias NCPolyElem NCPolyRingElem
@alias MPolyElem MPolyRingElem
@alias UnivPolyElem UniversalPolyRingElem
@alias LaurentPolyElem LaurentPolyRingElem
@alias LaurentMPolyElem LaurentMPolyRingElem

@alias AbsSeriesRing AbsPowerSeriesRing
@alias AbsSeriesElem AbsPowerSeriesRingElem
@alias RelSeriesRing RelPowerSeriesRing
@alias RelSeriesElem RelPowerSeriesRingElem
