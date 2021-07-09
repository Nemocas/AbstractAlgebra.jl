# Introduction

Maps in AbstractAlgebra model maps on sets $f : D \to C$ for some domain $D$
and codomain $C$, which have no real limitations except that elements of the
codomain and domain be represented by element objects in the system.

Maps $f : D \to C$ in AbstractAlgebra are modeled by Julia objects that are
able to be called on a single element $d \in D$ of the domain to yield an
element $f(d) \in C$ of the codomain. We say that the map is being applied.

Maps can be constructed from Julia functions, or they can be represented by
some other kind of data, e.g. a matrix, or built up from other maps.

Maps in AbstractAlgebra have a domain and codomain, can be applied, composed
and composed with the identity map (assuming its domain is compatible). Various
special kinds of map provide more functionality.

For example, there are functional maps which wrap a Julia function, cached
maps which cache values so they do not have to be recomputed each time they
are applied to the same inputs and various kinds of maps with inverses, e.g.
maps with sections, retractions and full inverses.

The map system uses a complex four parameter `Map` type, however various
helper functions are provided to make it easier to work with.

