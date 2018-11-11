# Constructing mathematical objects in AbstractAlgebra.jl

## Constructing objects in Julia

In Julia, one constructs objects of a given type by calling a type constructor. This is
simply a function with the same name as the type itself. For example, to construct a 
`BigInt` object in Julia, we simply call the `BigInt` constructor:

```julia
n = BigInt("1234567898765434567898765434567876543456787654567890")
```

Julia also uses constructors to convert between types. For example, to convert an `Int`
to a `BigInt`:

```julia
m = BigInt(123)
```

## How we construct objects in AbstractAlgebra.jl

As we explain in Appendix A, Julia types don't contain enough information to properly
model groups, rings, fields, etc. Instead of using types to construct objects, we use
special objects that we refer to as parent objects. They behave a lot like Julia types.

Consider the following simple example, to create a multiprecision integer:

```julia
n = ZZ("12345678765456787654567890987654567898765678909876567890")
```

Here `ZZ` is not a Julia type, but a callable object. However, for most purposes
one can think of such a parent object as though it were a type.

## Constructing parent objects

For more complicated groups, rings, fields, etc., one first needs to construct the
parent object before one can use it to construct element objects.

AbstractAlgebra.jl provides a set of functions for constructing such parent objects.
For example, to create a parent object for univariate polynomials over the integers,
we use the `PolynomialRing` parent object constructor.

```julia
R, x = PolynomialRing(ZZ, "x")
f = x^3 + 3x + 1
g = R(12)
```

In this example, `R` is the parent object and we use it to convert the `Int` value
$12$ to an element of the polynomial ring $\mathbb{Z}[x]$.

## List of parent object constructors

For convenience, we provide a list of all the parent object constructors in
AbstractAlgebra.jl and explain what mathematical domains they represent.

| Mathematics                      | AbstractAlgebra.jl constructor              |
|----------------------------------|---------------------------------------------|
| $R = \mathbb{Z}$                 | `R = ZZ`                                    |
| $R = \mathbb{Q}$                 | `R = QQ`                                    |
| $R = \mathbb{F}_{p}$             | `R = GF(p)`                                 |
| $R = \mathbb{Z}/n\mathbb{Z}$     | `R = ResidueRing(ZZ, n)`                    |
| $S = R[x]$                       | `S, x = PolynomialRing(R, "x")`             |
| $S = R[x, y]$                    | `S, (x, y) = PolynomialRing(R, ["x", "y"])` |
| $S = R[[x]]$ (to precision $n$)  | `S, x = PowerSeriesRing(R, n, "x")`         |
| $S = R((x))$ (to precision $n$)  | `S, x = LaurentSeriesRing(R, n, "x")`       |
| $S = K((x))$ (to precision $n$)  | `S, x = LaurentSeriesField(K, n, "x")`      |
| $S = \mbox{Frac}_R$              | `S = FractionField(R)`                      |
| $S = R/(f)$                      | `S = ResidueRing(R, f)`                     |
| $S = R/(f)$ (with $(f)$ maximal) | `S = ResidueField(R, f)`                     |
| $S = \mbox{Mat}_{m\times n}(R)$  | `S = MatrixSpace(R, m, n)`                  |
| $S = \mathbb{Q}[x]/(f)$          | `S, a = NumberField(f, "a")`                |

