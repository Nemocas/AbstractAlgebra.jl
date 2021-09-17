# Euclidean Ring Interface

If a ring provides a meaningful Euclidean structure such that a useful Euclidean
remainder can be computed practically, various additional functionality is provided
by AbstractAlgebra.jl for those rings. This functionality depends on the following
functions existing. An implementation must provide `divrem`, and the remaining
are optional as generic fallbacks exist.

```julia
divrem(f::MyElem, g::MyElem)
```

Return a pair `q, r` consisting of the Euclidean quotient and remainder of $f$
by $g$. A `DivideError` should be thrown if $g$ is zero.

```julia
mod(f::MyElem, g::MyElem)
```

Return the Euclidean remainder of $f$ by $g$. A `DivideError` should be thrown
if $g$ is zero.

```julia
div(f::MyElem, g::MyElem)
```

Return the Euclidean quotient of $f$ by $g$. A `DivideError` should be thrown
if $g$ is zero.

```julia
mulmod(f::MyElem, g::MyElem, m::MyElem)
```

Return $fg \pmod{m}$.

```julia
powermod(f::MyElem, e::Int, m::MyElem)
```

Return $f^e \pmod{m}$.

```julia
invmod(f::MyElem, m::MyElem)
```

Return the inverse of $f$ modulo $m$. If such an inverse doesn't exist, a
`NotInvertibleError` should be thrown.

```julia
divides(f::MyElem, g::MyElem)
```

Return a pair, `flag, q`, where `flag` is set to `true` if $g$ divides $f$, in which
case the quotient is set to the quotient or `flag` is set to `false` and the quotient
is set to zero in the same ring as $f$ and $g$.

```julia
remove(f::MyElem, p::MyElem)
```

Return a pair `v, q` where $p^v$ is the highest power of $p$ dividing $f$ and $q$ is
the cofactor after $f$ is divided by this power.

```julia
valuation(f::MyElem, p::MyElem)
```

Return `v` where $p^v$ is the highest power of $p$ dividing $f$.

```julia
gcd(f::MyElem, g::MyElem)
```

Return a greatest common divisor of $f$ and $g$. The return is expected to be
unit normalized such that if the return is a unit, that unit should be one.

```julia
lcm(f::MyElem, g::MyElem)
```

Return a least common multiple of $f$ and $g$.

```julia
gcdx(f::MyElem, g::MyElem)
```

Return a triple `d, s, t` such that $d = gcd(f, g)$ and $d = sf + tg$, with $s$
loosely reduced modulo $g/d$ and $t$ loosely reduced modulo $f/d$.

```julia
gcdinv(f::MyElem, g::MyElem)
```

Return a tuple `d, s` such that $d = gcd(f, g)$ and $s = (f/d)^{-1} \pmod{g/d}$. Note
that $d = 1$ iff $f$ is invertible modulo $g$, in which case $s = f^{-1} \pmod{g}$.
