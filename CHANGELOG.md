# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
tries to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.48.0] - 2025-12-XX
- TODO

## [0.47.0] - 2025-09-06

### BREAKING
- [#2132](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2132) Add `fraction_field_type` and define `fraction_field` method for fields
- [#2144](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2144) Remove `base_ring` methods returning `Union{}`
- [#2145](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2145) Add `identity_map` to excluded imports
- [#2150](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2150) Stop printing `1` in factorizations
- [#2153](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2153) Throw for `leading/trailing_coefficient` of the zero polynomial in multivariate setting

## [0.45.0] - 2025-05-15

### BREAKING
- [#2003](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2003) Remove `ignore_kwargs` option from `@attr`, that was deprecated in 0.44.4
- [#2002](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2002) Various changes to the conformance test setup:
  - Conformance tests are no longer loaded by including some file from the `test` directory, but are put into the submodule `ConformanceTests`, which only contains function stubs. The implementations get added by a package extension triggered by `Test`.
  - Thus, all code calling a function from the conformance tests (`test_mutating_op_*`, `test_*_interface(_recursive)`) needs to qualify the call with `ConformanceTests.`.
  - `test_elem` got renamed to `ConformanceTests.generate_element` and its implementations need to be moved form `test` to `src`.
- [#2006](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2006) Make `RationalFunctionFieldElem` immutable, remove a constructor method
- [#1809](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1809) Allow more inputs for `gens(::UniversalPolyRing, ...)` and change its return type from `Tuple` to `Vector`

### Added

- [#2072](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2072) Implement `is_perfect` for `RationalFunctionField`
- [#2081](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2081) Move `Infinity.jl` from Nemo
- [#2082](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2082) Move parts of Oscar's `MoveToAbstractAlgebra.jl` over

### Changed

- [#2071](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2071) Fix some cases where evaluating a multivariate polynomial returned an object of the wrong type


## [0.44.13] - 2025-05-07

### Changed

- [#2061](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2061) Fix a promote issue for `AbsMSeries` over Poly tower
- [#2065](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2065) Fix `swap_cols!` for `Matrix{T}`
- [#2068](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2068) Remove bad generic `sqrt` method that could lead to stack overflow
- [#2069](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2069) Fix `canonical_unit` methods returning zero

## [0.44.12] - 2025-04-25

### Added

- [#2039](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2039) Add generic `shift_left!` and `shift_right!` for univariate polynomials
- [#2047](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2047) Add `falling_factorial`

### Changed

- [#2058](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2058) Improve `to_univariate(::MPoly)`
- [#2056](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2056) Fix typo in `isomorphism(::Type{T}, ::T) where T <: Group` code
- [#2054](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2054) Fix stack overflow when constructing non-implemented views
- [#2048](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2048) Enable point views of matrices

## [0.44.11] - 2025-03-28

### Added

### Changed

- [#2042](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2042) Fix unwanted poly modification during evaluation

## [0.44.10] - 2025-03-19

### Added

- [#2035](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2035) Add stub for `is_equal_as_morphism`
- [#2036](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2036) Add stub for `isomorphism(::Type{<:Group}, ::Group)`

### Changed

- [#2037](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2037) Remove broken `is_unit` and `is_nilpotent` for `NCPolyRingElem`
- [#2033](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2033) Optimize `gen(::MPolyRing, ::Int)` for `:deglex` and `:degrevlex`

## [0.44.9] - 2025-03-10

### Added

- [#2020](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2020) Add `is_separable` for univariate polynomials
- [#2017](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2017) Add `is_positive_entry` and `is_negative_entry` for matrices
- [#2023](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2023) Add `reverse!` for polynomials
- [#2023](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2023) Add 2-argument `pow!` method delegating to the existing 3-argument version

### Changed

- [#2018](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2018) More intuitive polynomial evaluation
- [#2011](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2011) Add some missing `check_parent` in MPoly arithmetics
- [#2007](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2007) Wrap `MPoly` factorization for `UnivPoly`

## [0.44.8] - 2025-02-23

### Added

- [#2001](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2001) Add `minors_with_position` and `minors_iterator_with_position`
- [#2004](https://github.com/Nemocas/AbstractAlgebra.jl/pull/2004) Add `divides!`

### Changed

- [#1993](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1993) Allow creating `UniversalPolyRing` with variables
- [#1998](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1998) Fix `is_finite` for rational function fields

## [0.44.7] - 2025-02-17

### Added

### Changed

- [#1999](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1999) Delegate factorization for function fields to `FracField` over polynomial rings
- [#1999](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1999) Make `rand` for function fields produce reduced fractions

## [0.44.6] - 2025-02-13

### Added

- [#1991](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1991) Add `is_zero(::NCRing)` as an alias of `is_trivial`

### Changed

- [#1995](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1995) Add `is_finite` for rational function fields and allow omitting the variable name in `rational_function_field` (then `t` is used)
- [#1978](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1978) Import more `ConformanceTests` stuff into `TestExt`

## [0.44.5] - 2025-02-07

### Added

- [#1973](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1973) Enable multivariate factorization in characteristic 0 (at least over coefficient rings for which univariate factorization is implemented)
- [#1974](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1974) Add stub for `is_known`
- [#1962](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1962) Add changelog file
- [#1970](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1970) Add `content(::UnivPoly)`
- [#1954](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1954) Move conformance tests to a package extension `TestExt`

### Changed

- [#1982](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1982) Improve matrix documentation
- [#1983](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1983) Optimize `gens` for universal polynomials
- [#1977](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1977) Don't mention `charpoly_only` in docstrings, it is an internal helper
- [#1975](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1975) Improve documentation of `add_verbosity_scope`
- [#1970](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1970) Fix potentially wrong result in `content(::MatrixElem)`
- [#1970](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1970) Optimize a few `content` methods
- [#1971](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1971) Remove misleading `content` method
- [#1972](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1972) Fix docstring for `divides` promising too much

## [0.44.4] - 2025-01-22

### Added

### Changed

- [#1966](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1966) Ignore all kwargs in `@attr` (deprecating `ignore_kwargs`)
- [#1967](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1967) Fix the argument check in `remove` for `MPolyRingElem`


## [0.44.3] - 2025-01-17

### Added

- Add conformance tests for `MPolyRing`, and fix bugs in `^` and `is_unit` for `MPoly` [PR
  #1950](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1950)

- Add optional `ignore_kwargs` argument to `@attr` macro [PR
  #1958](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1958)

### Changed

- Don't require `base_ring` in ring conformance tests [issue
  #1944](https://github.com/Nemocas/AbstractAlgebra.jl/issues/1944) [PR
  #1946](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1946)

- Fix promotion for matrix-scalar operations [PR
  #1949](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1949)

- Some internal enhancements

## [0.44.2] - 2024-12-22

### Changed

- Make `is_exact_type` and `is_domain_type` more convenient by
  having them accept rings and ring elements (not just types) [PR
  #1942](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1942)

### Fixed

- Fix `is_unit` for univariate and multivariate polynomials over `Z/nZ` [issue
  #11](https://github.com/Nemocas/AbstractAlgebra.jl/issues/11) [PR
  #1933](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1933)

- Fix show method for MatSpace [PR #1934](https://github.com/Nemocas/AbstractAlgebra.jl/pull/1934)
