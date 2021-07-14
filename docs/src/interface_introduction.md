# Introduction

AbstractAlgebra defines a series of interfaces that can be extended with
new types that implement those interfaces. For example, if one were
implementing a new polynomial ring type, one would implement all of the
required functionality described in this chapter for the relevant
AbstractAlgebra interfaces. This would include the Ring Interface and the
Univariate Polynomial Ring Interface.

Once a new type implements all the required functionality, all the
corresponding generic functionality would then function automatically
for the new type.

One may then go on to implement some of the optional functionality for
performance if the provided generic functionality is insufficient.

AbstractAlgebra tries to provide all generic constructions recursively
so that one can have towers of generic constructions. This means that
new interfaces should generally only be added if they cooperate with all
the existing interfaces, at least so far as the theory exists to do so.

