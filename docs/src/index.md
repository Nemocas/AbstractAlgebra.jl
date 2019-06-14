# AbstractAlgebra.jl

## Introduction

AbstractAlgebra.jl is a computer algebra package for the Julia programming language, 
maintained by William Hart, Tommy Hofmann, Claus Fieker and Fredrik Johansson and
other interested contributors.

- [Source code](https://github.com/Nemocas/AbstractAlgebra.jl)
- [Online documentation](http://nemocas.github.io/AbstractAlgebra.jl)

AbstractAlgebra.jl grew out of the Nemo project after a number of requests from the
community for the pure Julia part of Nemo to be split off into a separate project. See
the Nemo website for more details about Nemo.
 
- [Nemo website](http://nemocas.org)

## Features

The features of AbstractAlgebra.jl include:

  - Use of Julia multiprecision integers and rationals
  - Finite fields (prime order, naive implementation only)
  - Number fields (naive implementation only)
  - Univariate polynomials
  - Multivariate polynomials
  - Relative and absolute power series
  - Laurent series
  - Fraction fields
  - Residue rings, including ``\mathbb{Z}/n\mathbb{Z}``
  - Matrices and linear algebra

All implementations are fully recursive and generic, so that one can build matrices
over polynomial rings, over a finite field, for example.

AbstractAlgebra.jl also provides a set of abstract types for Groups, Rings, Fields,
Modules and elements thereof, which allow external types to be made part of the
AbstractAlgebra.jl type hierarchy.

## Installation

To use AbstractAlgebra we require Julia 0.6 or higher. Please see
[http://julialang.org/downloads](http://julialang.org/downloads/) for instructions on 
how to obtain Julia for your system.

At the Julia prompt simply type

```
julia> Pkg.add("AbstractAlgebra")
```

## Quick start

Here are some examples of using AbstractAlgebra.jl.

This example makes use of multivariate polynomials.

```julia
using AbstractAlgebra

R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])

f = x + y + z + 1

p = f^20;

@time q = p*(p+1);
```

Here is an example using generic recursive ring constructions.

```julia
using AbstractAlgebra

R = GF(7)

S, y = PolynomialRing(R, "y")

T = ResidueRing(S, y^3 + 3y + 1)

U, z = PolynomialRing(T, "z")

f = (3y^2 + y + 2)*z^2 + (2*y^2 + 1)*z + 4y + 3;

g = (7y^2 - y + 7)*z^2 + (3y^2 + 1)*z + 2y + 1;

s = f^4;

t = (s + g)^4;

@time resultant(s, t)
```

Here is an example using matrices.

```julia
using AbstractAlgebra

R, x = PolynomialRing(ZZ, "x")

S = MatrixSpace(R, 10, 10)

M = rand(S, 0:3, -10:10);

@time det(M)
```

And here is an example with power series.

```julia
using AbstractAlgebra

R, x = QQ["x"]

S, t = PowerSeriesRing(R, 30, "t")

u = t + O(t^100)

@time divexact((u*exp(x*u)), (exp(u)-1));
```
