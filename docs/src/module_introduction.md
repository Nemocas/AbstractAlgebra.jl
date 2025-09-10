```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```
# Introduction

As with many generic constructions in AbstractAlgebra, the modules that are
provided in AbstractAlgebra itself work over a Euclidean domain. Moreover,
they are limited to finitely presented modules.

Free modules and vector spaces are provided over Euclidean domains and fields
respectively and then submodule, quotient module and direct sum module
constructions are possible recursively over these.

It's also possible to compute an invariant decomposition using the Smith
Normal Form.

The system also provides module homomorphisms and isomorphisms, building on
top of the map interface.

As for rings and fields, modules follow an interface which other modules are
expected to follow. However, very little generic functionality is provided
automatically once this interface is implemented by a new module type.

The purpose of the module interface is simply to encourage uniformity in
the module interfaces of systems that build on AbstractAlgebra. Of course
modules are so diverse that this is a very loosely defined interface to
accommodate the diversity of possible representations and implementations.




