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

In fact, the definition of the number field constructor is currently given in
`src/generic/ResidueField.jl` and no type is defined for a number field. The definition
mainly exists for testing purposes. It may later be replaced by a more standard
implementation. For a more fully fleshed out number field implementation (based on a
very high performance C library), see Nemo.jl.

## Number field constructors

In order to construct number fields in AbstractAlgebra.jl, one must first construct the
field itself. This is accomplished with the following constructor.

```julia
NumberField(f::AbstractAlgebra.Generic.Poly{Rational{BigInt}}, s::AbstractString, t = "\$"; cached = true)
```

Given an irreducible defining polynomial $f$ in $\mathbb{Q}[x]$, return a tuple $(K, x)$
consisting of the number field defined by that polynomial and a generator. The string
fields are currently ignored, but are reserved for future use.

Currently the generator of the number field prints the same way as the variable in
$\mathbb{Q}[x]$.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> K, a = NumberField(x^3 + 3x + 1, "a")
(Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1, x)

julia> f = a^2 + 2a + 7
x^2+2//1*x+7//1

```

## Basic field functionality

The number field module in AbstractAlgebra.jl implements the full Field and ResidueRing
interfaces.

