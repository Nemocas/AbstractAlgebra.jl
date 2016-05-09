# Nemo

## Introduction

Nemo is a computer algebra package for the Julia programming language, maintained by William Hart with code by William Hart, Tommy Hofmann, Claus Fieker, Fredrik Johansson, Oleksandr Motsak and other contributors.

- [https://github.com/Nemocas/Nemo.jl](https://github.com/Nemocas/Nemo.jl) (Source code)
- [http://nemo.readthedocs.org/en/latest/](http://nemo.readthedocs.org/en/latest/) (Online documentation)

The features of Nemo so far include:

  - Multiprecision integers and rationals
  - Integers modulo n
  - p-adic numbers
  - Finite fields (prime and non-prime order)
  - Number field arithmetic
  - Maximal orders of number fields
  - Arithmetic of ideals in maximal orders
  - Arbitrary precision real and complex balls
  - Generic polynomials, power series, fraction fields, residue rings and matrices

## Installation

To use Nemo we require Julia 0.4 or higher. Please see [http://julialang.org/downloads](http://julialang.org/downloads/) for instructions on how to obtain julia for your system.

At the Julia prompt simply type

```
julia> Pkg.add("Nemo")
julia> Pkg.build("Nemo")
```


