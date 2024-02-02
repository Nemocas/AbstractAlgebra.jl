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

The attributes system can be utilized to change the way certain objects are printed.
We provide macros `@show_special` and `@show_name` for this purpose, both are
called with the same argument
as `show`: an `IO`-object and the object itself. Both are supposed to be
used within the usual `show` function:
```
function show(io::IO, A::MyObj)
   @show_name(io, A)
   @show_special(io, A)

   ... usual stuff
```  

`@show_special` checks if an attribute `:show` is present. If so, it has to be
a function taking `IO` and the object. This is then called instead of the usual
`show` function.

`@show_name` will check if there is a variable in global (`Main` module) namespace
with value bound to the object. In compact printing mode, the name is then shown
instead of the object.

Note: if the object is stored in several variable, the first one will be used. Also
the name, once used for printing, is stored in the object - hence will not change
anymore.

## Advanced printing

To facilitate printing of nested mathematical structures, we provide a modified
`IOCustom` object, that supports indentation and decapitalization.

### Example

We illustrate this with an example

```
struct A{T}
  x::T
end

function Base.show(io::IO, a::A)
  io = AbstractAlgebra.pretty(io)
  println(io, "Something of type A")
  print(io, AbstractAlgebra.Indent(), "over ", AbstractAlgebra.Lowercase(), a.x)
  print(io, AbstractAlgebra.Dedent()) # don't forget to undo the indentation!
end

struct B
end

function Base.show(io::IO, b::B)
  io = AbstractAlgebra.pretty(io)
  print(io, LowercaseOff(), "Hilbert thing")
end
```

At the REPL, this will then be printed as follows:
```
julia> A(2)
Something of type A
  over 2

julia> A(A(2))
Something of type A
  over something of type A
    over 2

julia> A(B())
Something of type A
  over Hilbert thing
```

### Documentation

```@docs
AbstractAlgebra.pretty
AbstractAlgebra.Indent
AbstractAlgebra.Dedent
AbstractAlgebra.Lowercase
AbstractAlgebra.LowercaseOff
```
