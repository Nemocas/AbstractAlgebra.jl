# Noncommutative ring Interface

AbstractAlgebra.jl supports commutative rings through its `Ring` interface.
In this section we describe the corresponding interface for noncommutative
rings. The two interfaces are very similar in terms of required functionality,
and so we mainly document the differences here.

Noncommutative rings can be supported through the abstract types
`NCRing` and `NCRingElem`. Note that we have
`Ring <: NCRing`, etc., so the interface here
should more correctly be called the Not-necessarily-Commutative-ring interface.

However, the fact remains that if one wishes to implement a noncommutative
ring, one should make its type belong to `NCRing` but not to
`Ring`. Therefore it is not too much of a mistake to think of
the `NCRing` interface as being for noncommutative rings.

## Types

As for the Ring interface, most noncommutative rings must supply two types:
  - a type for the parent object (representing the ring itself)
  - a type for elements of that ring

The parent type must belong to `NCRing` and the element type
must belong to `NCRingElem`. Of course, the types may belong
to these abstract types transitively via an intermediate abstract type.

Also as for the Ring interface, it is advised to make the types of generic
parameterised rings that belong to `NCRing` and `NCRingElem` depend on the type
of the elements of that parameter ring.

## NCRingElement type union

As for the Ring interface, the NCRing interface provides a union type
`NCRingElement` in `src/julia/JuliaTypes.jl` which is a union
of `NCRingElem` and the Julia types `Integer`, `Rational`
and `AbstractFloat`.
                                                                                                    
Most of the generic code in AbstractAlgebra for general rings makes use of the
union type `NCRingElement` instead of `NCRingElem`
so that the generic functions also accept the Julia Base ring types.
              
As per usual, one may need to implement one ad hoc binary operation for each
concrete type belonging to `NCRingElement` to avoid ambiguity
warnings.

## Parent object caches

Parent object caches for the NCRing interface operate as per the Ring interface.

## Required functions for all rings

Generic functions may only rely on required functionality for the NCRing
interface, which must be implemented by all noncommutative rings.

Most of this required functionality is the same as for the Ring interface, so
we refer the reader there for details, with the following modifications.

We give this interface for fictitious types `MyParent` for the type of the ring parent
object `R` and `MyElem` for the type of the elements of the ring.

### Exact division

```julia
divexact_left(f::MyElem, g::MyElem)
divexact_right(f::MyElem, g::MyElem)
```

If $f = ga$ for some $a$ in the ring, the function `divexact_left(f, g)` returns `a`. If
$f = ag$ then `divexact_right(f, g)` returns `a`. A `DivideError()` should be thrown
if division is by zero. If no exact quotient exists or an impossible inverse is
unavoidably encountered, an error should be thrown.

