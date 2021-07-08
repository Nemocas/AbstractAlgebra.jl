# Introduction

A rich ring hierarchy is provided, supporting both commutative and
noncommutative rings.

A number of basic rings are provided, such as the integers, integers
mod `n` and numerous fields.

A recursive rings implementation is then built on top of the basic
rings via a number of generic ring constructions. These include
univariate and multivariate polynomials and power series, univariate
Laurent and Puiseux series, residue rings, matrix algebras, etc.

Where possible, these constructions can be built on top of one another
in generic towers.

The ring hierarchy can be extended by implementing new rings to follow
one or more ring interfaces. Generic functionality provided by the
system is then automatically available for the new rings. These
implementations can either be generic or can be specialised
implementations provided by, for example, a C library.

In most cases, the interfaces consist of a set of constructors and
functions that must be implemented to satisfy the interface. These are
the functions that the generic code relies on being available.


