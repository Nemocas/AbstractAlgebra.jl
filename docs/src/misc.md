```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
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
is_attribute_storing
is_attribute_storing_type
```

## Advanced printing

### Self-given names

We provide macros `@show_name`, `@show_special` and `@show_special_elem` to 
change the way certain objects are printed.

In compact and terse printing mode, `@show_name` tries to determine
a suitable name to print instead of the object (see [`AbstractAlgebra.get_name`](@ref)).

`@show_special` checks if an attribute `:show` is present. If so, it has to be
a function taking `IO`, optionally a MIME-type, and the object.
This is then called instead of the usual `show` function.

Similarly, `@show_special_elem` checks if an attribute `:show_elem` is present in the object's
parent. The semantics are the same as for `@show_special`.

All are supposed to be used within the usual `show` function, where `@show_special_elem`
is only relevant for element types of algebraic structures.
```julia
@attributes MyObj

function show(io::IO, A::MyObj)
   @show_name(io, A)
   @show_special(io, A)

   # ... usual stuff
end

function show(io::IO, mime::MIME"text/plain", A::MyObj)
   @show_name(io, A)
   @show_special(io, mime, A)

   # ... usual stuff
end

function show(io::IO, A::MyObjElem)
   @show_name(io, A)
   @show_special_elem(io, A)

   # ... usual stuff
end

function show(io::IO, mime::MIME"text/plain", A::MyObjElem)
   @show_name(io, A)
   @show_special_elem(io, mime, A)

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
AbstractAlgebra.PrettyPrinting.find_name
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

```@docs
AbstractAlgebra.terse
AbstractAlgebra.is_terse
```

## Linear solving interface for developers

AbstractAlgebra has a generic interface for linear solving and we describe here
how one may extend this interface. For the user-facing functionality of linear
solving, see [Linear Solving](@ref solving_chapter).

Notice that the functionality is implemented in the module
`AbstractAlgebra.Solve` and the internal functions are not exported from there.

### Matrix normal forms

To distinguish between different algorithms, we use type traits of abstract type
`MatrixNormalFormTrait` which usually correspond to a certain matrix normal
form.
The available algorithms/normal forms are
* `HowellFormTrait`: uses a Howell form;
* `HermiteFormTrait`: uses a Hermite normal form;
* `RREFTrait`: uses a row-reduced echelon form over fields;
* `LUTrait`: uses a LU factoring of the matrix;
* `FFLUTrait`: uses a "fraction-free" LU factoring of the matrix over fraction
    fields;
* `MatrixInterpolateTrait`: uses interpolation of polynomials for fraction
    fields of polynomial rings.

To select a normal form type for rings of type `NewRing`, implement the function
```julia
Solve.matrix_normal_form_type(::NewRing) = Bla()
```
where `Bla <: MatrixNormalFormTrait`.
A new type trait can be added via
```julia
struct NewTrait <: Solve.MatrixNormalFormTrait end
```

### Internal solving functionality

If a new ring type `NewRing` can make use of one of the available
`MatrixNormalFormTrait`s, then it suffices to specify this normal form as
described above to use the generic solving functionality. (However, for example
`HermiteFormTrait` requires that the function
`hermite_form_with_transformation` is implemented.)

For a new trait `NewTrait <: MatrixNormalFormTrait`, one needs to implement the
function
```julia
Solve._can_solve_internal_no_check(
  ::NewTrait, A::MatElem{T}, b::MatElem{T}, task::Symbol; side::Symbol = :left
  ) where T
```
Inside this function, one can assume that `A` and `b` have the same base ring
and have compatible dimensions. Further, `task` and `side` are set to "legal"
options. (All this is checked in `Solve._can_solve_internal`.)
This function should then (try to) solve `Ax = b` (`side == :right`) or `xA = b`
(`side == :left`) possibly with kernel.
The function must always return a tuple `(::Bool, ::MatElem{T}, ::MatElem{T})`
consisting of:
* `true`/`false` whether a solution exists or not
* the solution (or a placeholder if no solution exists or a solution is not requested)
* the kernel (or a placeholder if the kernel is not requested)

The input `task` may be:
* `:only_check`: Only test whether there is a solution, the second
  and third return value are only for type stability;
* `:with_solution`: Compute a solution, if it exists, the last return value is
  only for type stability;
* `:with_kernel`: Compute a solution and a kernel.

One should further implement the function
```julia
kernel(::NewTrait, A::MatElem; side::Symbol = :left)
```
which computes a left (or right) kernel of `A`.

### Internal solve context functionality

To efficiently solve several linear systems with the same matrix `A`, we provide
the "solve contexts objects" of type `Solve.SolveCtx`.
These can be extended for a ring of type `NewRing` as follows.

#### Solve context type

For a new ring type, one may have to define the type parameters of a `Solve.SolveCtx`
object.
First of all, one needs to implement the function

```julia
function Solve.solve_context_type(::NewRing)
  return Solve.solve_context_type(::NormalFormTrait, elem_type(NewRing))
end
```

to pick a `MatrixNormalFormTrait`.

Usually, nothing else should be necessary. However, if for example the normal
form of a matrix does not live over the same ring as the matrix itself, one
might also need to implement

```julia
function Solve.solve_context_type(NF::NormalFormTrait, T::Type{NewRingElem})
  return Solve.SolveCtx{T, typeof(NF), MatType, RedMatType, TranspMatType}
end
```

where `MatType` is the dense matrix type over `NewRing`, `RedMatType` the type
of a matrix in reduced/normal form and `TranspMatType` the type of the
reduced/normal form of the transposed matrix.

#### Initialization

To initialize the solve context functionality for a new normal form `NewTrait`,
one needs to implement the functions

```julia
Solve._init_reduce(C::Solve.SolveCtx{T, NewTrait}) where T
Solve._init_reduce_transpose(C::Solve.SolveCtx{T, NewTrait}) where T
```

These should fill the corresponding fields of the solve context `C` with a
"reduced matrix" (that is, a matrix in normal form) of `matrix(C)`, respectively
`transpose(matrix(C))`, and other information necessary to solve a linear system.
The fields can be accessed via `reduced_matrix`, `reduced_matrix_of_transpose`,
etc. New fields may also be added via attributes.

#### Internal solving functionality

As above, one finally needs to implement the functions

```julia
Solve._can_solve_internal_no_check(
  ::NewTrait, C::Solve.SolveCtx{T, NewTrait}, b::MatElem{T}, task::Symbol;
  side::Symbol = :left
  ) where T
```

and

```julia
kernel(::NewTrait, C::Solve.SolveCtx{T, NewTrait}; side::Symbol = :left)
```
