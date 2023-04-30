```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Number fields

AbstractAlgebra.jl provides a very naive implementation of number fields. This allows
arithmetic in algebraic number fields, which are currently modeled as $\mathbb{Q}[x]$
modulo an irreducible polynomial, i.e. as a residue field.

The definition of the number field constructor is given in
`src/generic/NumberField.jl` but no type is defined for a number field. The definition
mainly exists for testing purposes. It may later be replaced by a more standard
implementation. For a more fully fleshed out number field implementation (based on a
very high performance C library), see Nemo.jl.

## Number field constructors

In order to construct number fields in AbstractAlgebra.jl, one must first construct the
field itself. This is accomplished with the following constructor.

```julia
number_field(f::Generic.Poly{Rational{BigInt}}, s::VarName, t = "\$"; cached::Bool = true)
```

Given an irreducible defining polynomial $f$ in $\mathbb{Q}[x]$, return a tuple $(K, x)$
consisting of the number field defined by that polynomial and a generator. The string
fields are currently ignored, but are reserved for future use.

Currently the generator of the number field prints the same way as the variable in
$\mathbb{Q}[x]$.

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> K, a = number_field(x^3 + 3x + 1, "a")
(Residue field of univariate polynomial ring modulo x^3 + 3*x + 1, x)

julia> f = a^2 + 2a + 7
x^2 + 2*x + 7

```

## Basic field functionality

The number field module in AbstractAlgebra.jl implements the full Field and residue_ring
interfaces.

