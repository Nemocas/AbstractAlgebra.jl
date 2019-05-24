# AbstractAlgebra

[![Build Status](https://travis-ci.org/Nemocas/AbstractAlgebra.jl.svg?branch=master)](https://travis-ci.org/Nemocas/AbstractAlgebra.jl) [![Build status](https://ci.appveyor.com/api/projects/status/1w9ninmoidxkxshp/branch/master?svg=true)](https://ci.appveyor.com/project/thofma/abstractalgebra-jl/branch/master)

AbstractAlgebra is a pure Julia package for computational abstract algebra. It grew out of the Nemo project and provides all of the abstract types and generic implementations that Nemo relies on.

It is currently developed by William Hart, Tommy Hofmann, Fredrik Johansson,
Claus Fieker with contributions from others.

AbstractAlgebra currently provides:

* Generic polynomial rings, matrix spaces, fraction fields, residue rings, relative and absolute power series, Laurent series
* Finite fields, integers, rationals, permutations and characters, number fields

Documentation can be found at the following link:

* [https://nemocas.github.io/AbstractAlgebra.jl/latest/index.html](https://nemocas.github.io/AbstractAlgebra.jl/latest/index.html)

Projects that depend on AbstractAlgebra include:

* Nemo.jl [http://nemocas.org](http://nemocas.org/) (optimised implementations of specific rings provided by the Flint, Arb and Antic C libraries)
* Hecke.jl [https://github.com/thofma/Hecke.jl](https://github.com/thofma/Hecke.jl) (algebraic number theory)
* Singular.jl [https://github.com/oscar-system/Singular.jl](https://github.com/oscar-system/Singular.jl) (polynomial rings and ideals, Groebner bases and computer algebra provided by the Singular C++ library)

