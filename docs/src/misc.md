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

## Advanced printing

### Self-given names

We provide macros `@show_name`, `@show_special` and `@show_special_elem` to 
change the way certain objects are printed.

In compact and supercompact printing mode, `@show_name` tries to determine
a suitable name to print instead of the object (see [`AbstractAlgebra.get_name`](@ref)).

`@show_special` checks if an attribute `:show` is present. If so, it has to be
a function taking `IO`, optionally a MIME-type, and the object.
This is then called instead of the usual `show` function.

Similarly, `@show_special_elem` checks if an attribute `:show_elem` is present in the object's
parent. The semantics are the same as for `@show_special`.

All are supposed to be used within the usual `show` function, where `@show_special_elem`
is only relevant for element types of algebraic structures.
```julia
function show(io::IO, A::MyObj)
   @show_name(io, A)
   @show_special(io, A)

   # ... usual stuff
end

function show(io::IO, ::MIME"text/plain", A::MyObj)
   @show_name(io, A)
   @show_special(io, MIME"text/plain"(), A)

   # ... usual stuff
end
```

#### Documentation

```@docs
AbstractAlgebra.@show_special
AbstractAlgebra.@show_special_elem
AbstractAlgebra.@show_name
AbstractAlgebra.get_name
AbstractAlgebra.set_name!
AbstractAlgebra.extra_name
AbstractAlgebra.find_name
```

### Indentation and Decapitalization

To facilitate printing of nested mathematical structures, we provide a modified
`IOCustom` object, that supports indentation and decapitalization.

#### Example

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

#### Documentation

```@docs
AbstractAlgebra.pretty
AbstractAlgebra.Indent
AbstractAlgebra.Dedent
AbstractAlgebra.Lowercase
AbstractAlgebra.LowercaseOff
```
