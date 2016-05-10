# Nemo

## What is Nemo?

Nemo is a computer algebra package for the Julia programming language. Our aim is to provide a highly
performant computer algebra package covering

  - Commutative Algebra
  - Number Theory
  - Group Theory

Nemo consists of wrappers of specialised C/C++ libraries:

  - Flint    [http://flintlib.org/]
  - Arb      [http://fredrikj.net/arb/]
  - Antic    [https://github.com/wbhart/antic/]
  - Singular [https://www.singular.uni-kl.de/]
  - Pari     [http://pari.math.u-bordeaux.fr/]

It will also eventually provide interfaces to interpreted library code from other computer algebra
systems such as Gap and Singular.

Nemo also provides implementations of generic algorithms and mathematical data structures. So far the
fully recursive constructions include

  - Univariate polynomial rings
  - Power series rings
  - Residue rings (modulo principal ideals)
  - Fraction fields
  - Matrices
