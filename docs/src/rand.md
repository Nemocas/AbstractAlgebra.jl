```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Random interface

AbstractAlgebra makes use of the Julia
[Random interface](https://docs.julialang.org/en/v1/stdlib/Random/) for random
generation.

In addition we make use of an experimental package
[RandomExtensions.jl](https://github.com/JuliaRandom/RandomExtensions.jl) for
extending the random interface in Julia.

The latter is required because some of our types require more than one
argument to specify how to randomise them.

The usual way of generating random values that Julia and these extensions
provide would look as follows:

```julia
julia> using AbstractAlgebra

julia> using Random

julia> using RandomExtensions

julia> S, x = polynomial_ring(ZZ, :x)
(Univariate Polynomial Ring in x over Integers, x)

julia> rand(Random.GLOBAL_RNG, make(S, 1:3, -10:10))
-5*x + 4
```

This example generates a polynomial over the integers with degree in the range
1 to 3 and with coefficients in the range -10 to 10.

In addition we implement shortened versions for ease of use which don't
require creating a make instance or passing in the standard RNG.

```julia
julia> using AbstractAlgebra

julia> S, x = polynomial_ring(ZZ, :x)
(Univariate Polynomial Ring in x over Integers, x)

julia> rand(S, 1:3, -10:10)
-5*x + 4
```

Because rings can be constructed over other rings in a tower, all of this is
supported by defining `RandomExtensions.make` instances that break the various
levels of the ring down into separate make instances.

For example, here is the implementation of `make` for polynomial rings such as
the above:

```julia
function RandomExtensions.make(S::PolyRing, deg_range::AbstractUnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, deg_range, vs[1]) # forward to default Make constructor
   else
      Make(S, deg_range, make(R, vs...))
   end
end
```

As you can see, it has two cases. The first is where this invocation of make is
already at the bottom of the tower of rings, in which case it just forwards to
the default `Make` constructor.

The second case expects that we are higher up in the tower of rings and that
`make` needs to be broken up (recursively) into the part that deals with the
ring level we are at and the level that deals with the base ring.

To help `make` we tell it the type of object we are hoping to randomly
generate.

```julia
RandomExtensions.maketype(S::PolyRing, dr::AbstractUnitRange{Int}, _) = elem_type(S)
```

Finally we implement the actual random generation itself.

```julia
# define rand for make(S, deg_range, v)
function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make3{<:RingElement, <:PolyRing, <:AbstractUnitRange{Int}}})
   S, deg_range, v = sp[][1:end]
   R = base_ring(S)
   f = S()
   x = gen(S)
   # degree -1 is zero polynomial
   deg = rand(rng, deg_range)
   if deg == -1
      return f
   end
   for i = 0:deg - 1
      f += rand(rng, v)*x^i
   end
   # ensure leading coefficient is nonzero
   c = R()
   while iszero(c)
      c = rand(rng, v)
   end
   f += c*x^deg
   return f
end
```

Note that when generating random elements of the base ring for example, one
should use the random number generator `rng` that is passed in.

As mentioned above, we define a simplified random generator that saves the user
having to create make instances.

```julia
rand(rng::AbstractRNG, S::PolyRing, deg_range::AbstractUnitRange{Int}, v...) =
   rand(rng, make(S, deg_range, v...))

rand(S::PolyRing, degs, v...) = rand(Random.GLOBAL_RNG, S, degs, v...)
```

To test whether a random generator is working properly, the `test_rand` function
exists in the AbstractAlgebra test submodule in the file `test/runtests.jl`.
For example, in AbstractAlgebra test code:

```julia
using Test

R, x = polynomial_ring(ZZ, :x)

test_rand(R, -1:10, -10:10)
```

In general, we try to use `UnitRange`'s to specify how 'big' we want the
random instance to be, e.g. the range of degrees a polynomial could take,
the range random integers could lie in, etc. The objective is to make it
easy for the user to control the 'size' of random values in test code.

