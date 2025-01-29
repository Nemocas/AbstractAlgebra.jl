# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
tries to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Added

- Started keeping a changelog!

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
