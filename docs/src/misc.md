```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Miscellaneous

## Printing options

AbstractAlgebra supports printing to LaTeX using the MIME type "text/latex". To
enable LaTeX rendering in Jupyter notebooks and query for the current state,
use the following functions:

```@docs
set_html_as_latex
get_html_as_latex
```

## Updating the type diagrams

Updating the diagrams of the documentation can be done by modifying and running
the script `docs/create_type_diagrams.jl`. Note that this requires the package `Kroki`.


## Attributes

Often it is desirable to have a flexible way to attach additional data to
mathematical structures such as groups, rings, fields, etc. beyond what the
original implementation covers. To facilitate this, we provide an *attributes*
system: for objects of suitable types, one may use `set_attribute!` to attach
key-value pairs to the object, and query them using `has_attribute`,
`get_attribute` and `get_attribute!`.

Attributes are supported for all singletons (i.e., instances of an empty
`struct` type), as well as for instances of mutable struct type for which
attribute storage was enabled. There are two ways to enable attribute storage
for such types:

1. By applying `@attributes` to a mutable struct declaration, storage is
   reserved inside that struct type itself (this increases the size of each
   struct by 8 bytes if no attributes are set).
2. By applying `@attributes` to the name of a mutable struct type, methods are
   installed which store attributes to instances of the type in a
   `WeakKeyDict` outside the struct.

```@docs
@attributes
@attr
has_attribute
get_attribute
get_attribute!
set_attribute!
```
