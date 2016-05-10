# Types in Nemo

## Introduction

Julia provides two levels of types that we make use of

  - abstract types
  - concrete types

Concrete types are just like the usual types everyone is familiar with from C or C++.

Abstract types can be thought of as collections of types. They are used when writing generic functions
that should work for any type in the given collection.

To write a generic function that accepts any type in a given collection of types, we first create an
abstract type. Then we create the individual concrete types that belong to that abstract type. A generic
function is then constructed with a type parameter, `T` say, similar to a template class in C++. The main
difference is that we can specify which abstract type our type parameter `T` must belong to.

We use the symbol <: in Julia to determine that a given type belongs to a given abstract type. 

Here is some Julia code illustrating this. We create an abstract type called `Shape` and two user defined
concrete types `square` and `circle` belonging to `Shape`. We then show how to write methods for the
concrete types and also how to write a generic function for any type `T` belonging to the abstract type
`Shape`.

```
abstract Shape

type square <: Shape
   width::Int
   border_thickness::Int
end

type circle <: Shape
   centre::Tuple{Int, Int}
   radius::Int
   border_thickness::Int
end

function area(s::square)
   return s.width^2
end

function area(s::circle)
   return pi*s.radius^2
end

function border_thickness{T <: Shape}(s::T)
   return s.border_thickness
end

s = square(3, 1)
c = circle((3, 4), 2, 2)

area(s)
area(c)
border_thickness(s)
```

## The abstract type hierarchy in Nemo

Abstract types in Julia can also belong to one another in a hierarchy. For example, the `Nemo.Field`
abstract type belongs to the `Nemo.Ring` abstract type. This means that any object representing a field
in Nemo has a type belonging to both `Nemo.Field` and `Nemo.Ring`, the latter automatically because of
the inclusion `Nemo.Field <: Nemo.Ring`.

In Nemo we also distinguish between the elements of a field, say, and the field itself. For example,
we have an object of type `PolynomialRing` to model a generic polynomial ring, and elements of that 
polynomial ring would have type `Poly`. 

In order to model this distinction between elements and the domains they belong to, Nemo has two main
branches in its abstract type hierarchy, as shown in the following diagram.

![alt text](/types.png "Abstract type hierarchy")

All objects in Nemo, whether they represent rings, fields, groups, sets on the one hand, or ring
elements, field elements, etc. on the other hand, have concrete types that belong to one of the abstract
types shown above.

