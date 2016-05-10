# Types in Nemo

## Introduction

Julia provides two levels of types that we make use of

  - abstract types
  - concrete types

Concrete types are just like the usual types everyone is familiar with from C or C++.

Abstract types can be thought of as collections of types. They are used when writing generic functions
that should work for any type in the given collection.

To write a generic function that accepts any type in a given collection of types, we create a type
parameter, similar to the way template functions are created in C++, and specify that the type
parameter belongs to the given abstract type. We use the symbol <: in Julia to determine if a given
type belongs to a given abstract type. 

Here is some Julia code illustrating this. We create an abstract type called Shape and two user defined
concrete types square and circle. We then show how to write methods for various combinations of these.

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
   