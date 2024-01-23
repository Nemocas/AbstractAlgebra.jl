# for backwards compatibility
Base.@deprecate_binding iscompatible is_compatible
Base.@deprecate_binding isconstant is_constant
Base.@deprecate_binding isdegree is_degree
Base.@deprecate_binding isdivisible_by is_divisible_by
Base.@deprecate_binding isdomain_type is_domain_type
Base.@deprecate_binding isexact_type is_exact_type
Base.@deprecate_binding isgen is_gen
Base.@deprecate_binding ishessenberg is_hessenberg
Base.@deprecate_binding ishnf is_hnf
Base.@deprecate_binding ishomogeneous is_homogeneous
Base.@deprecate_binding isinvertible is_invertible
Base.@deprecate_binding isinvertible_with_inverse is_invertible_with_inverse
Base.@deprecate_binding isirreducible is_irreducible
Base.@deprecate_binding isisomorphic is_isomorphic
Base.@deprecate_binding ismonic is_monic
Base.@deprecate_binding ismonomial is_monomial
Base.@deprecate_binding ismonomial_recursive is_monomial_recursive
Base.@deprecate_binding isnegative is_negative
Base.@deprecate_binding ispopov is_popov
Base.@deprecate_binding ispower is_power
Base.@deprecate_binding isprobable_prime is_probable_prime
Base.@deprecate_binding isreverse is_reverse
Base.@deprecate_binding isrimhook is_rimhook
Base.@deprecate_binding isrref is_rref
Base.@deprecate_binding issnf is_snf
Base.@deprecate_binding issquare is_square
Base.@deprecate_binding issquare_with_sqrt is_square_with_sqrt
Base.@deprecate_binding issquarefree is_squarefree
Base.@deprecate_binding issubmodule is_submodule
Base.@deprecate_binding isterm is_term
Base.@deprecate_binding isterm_recursive is_term_recursive
Base.@deprecate_binding isunit is_unit
Base.@deprecate_binding isunivariate is_univariate
Base.@deprecate_binding isweak_popov is_weak_popov
Base.@deprecate_binding iszero_column is_zero_column
Base.@deprecate_binding iszero_row is_zero_row
Base.@deprecate_binding Localization localization

#
# Some deprecation aliases follow. They use Base.@deprecate_binding instead of @deprecate
# because there are still packages adding methods for them, which would break
# if we used @deprecate. Once those packages are updated, we can switch to
# using @deprecate
#

# Deprecated in 0.27.*

Base.@deprecate_binding MatrixSpace matrix_space

Base.@deprecate_binding PolynomialRing polynomial_ring
Base.@deprecate_binding PowerSeriesRing power_series_ring
Base.@deprecate_binding LaurentPolynomialRing laurent_polynomial_ring
Base.@deprecate_binding LaurentSeriesRing laurent_series_ring
Base.@deprecate_binding LaurentSeriesField laurent_series_field
Base.@deprecate_binding FractionField fraction_field

Base.@deprecate_binding FunctionField function_field
Base.@deprecate_binding RationalFunctionField rational_function_field

Base.@deprecate_binding PolyElem PolyRingElem
Base.@deprecate_binding NCPolyElem NCPolyRingElem
Base.@deprecate_binding MPolyElem MPolyRingElem
Base.@deprecate_binding UnivPolyElem UniversalPolyRingElem
Base.@deprecate_binding LaurentPolyElem LaurentPolyRingElem
Base.@deprecate_binding LaurentMPolyElem LaurentMPolyRingElem

Base.@deprecate_binding UnivPolyRing UniversalPolyRing

Base.@deprecate_binding AbsSeriesRing AbsPowerSeriesRing
Base.@deprecate_binding AbsSeriesElem AbsPowerSeriesRingElem
Base.@deprecate_binding RelSeriesRing RelPowerSeriesRing
Base.@deprecate_binding RelSeriesElem RelPowerSeriesRingElem

Base.@deprecate_binding FreeAssociativeAlgebra free_associative_algebra

Base.@deprecate_binding TotalRingOfFractions total_ring_of_fractions

Base.@deprecate_binding DirectSum direct_sum

# Deprecated in 0.34.*
Base.@deprecate_binding PuiseuxSeriesField puiseux_series_field
Base.@deprecate_binding PuiseuxSeriesRing puiseux_series_ring
Base.@deprecate_binding istrivial is_trivial
