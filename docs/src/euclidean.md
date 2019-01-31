# Euclidean Ring Interface

If a ring provides a meaningful Euclidean structure such that a useful Euclidean
remainder can be computed practically, various additional functionality is provided
by AbstractAlgebra.jl for those rings. This functionality depends on the following
functions existing.

```julia
mod(f::MyElem, g::MyElem)
```

Returns the Euclidean remainder of $f$ by $g$. A `DivideError()` should be thrown if
$g$ is zero. An error should be thrown if an impossible inverse is encountered.

```julia
divrem(f::MyElem, g::MyElem)
```

Returns a pair `q, r` consisting of the Euclidean quotient and remainder of $f$ by $g$.
A `DivideError` should be thrown if $g$ is zero. An error should be thrown if an
impossible inverse is encountered.

```julia
div(f::MyElem, g::MyElem)
```

Returns the Euclidean quotient of $f$ by $g$. A `DivideError` should be thrown if $g$
is zero. An error should be thrown if an impossible inverse is encountered.

```julia
mulmod(f::MyElem, g::MyElem, m::MyElem)
```

Returns $fg \pmod{m}$.

```julia
powmod(f::MyElem, e::Int, m::MyElem)
```

Returns $f^e \pmod{m}$.

```julia
invmod(f::MyElem, m::MyElem)
```

Returns the inverse of $f$ modulo $m$. If such an inverse doesn't exist, an impossible
inverse error should be thrown.

```julia
divides(f::MyElem, g::MyElem)
```

Returns a pair, `flag, q`, where `flag` is set to `true` if $g$ divides $f$, in which
case the quotient is set to the quotient or `flag` is set to `false` and the quotient
is set to zero in the same ring as $f$ and $g$.

```julia
remove(f::MyElem, p::MyElem)
```

Returns a pair `v, q` where $p^v$ is the highest power of $p$ dividing $f$ and $q$ is
the cofactor after $f$ is divided by this power.

```julia
valuation(f::MyElem, p::MyElem)
```

Returns `v` where $p^v$ is the highest power of $p$ dividing $f$.

```julia
gcd(f::MyElem, g::MyElem)
```

Returns a greatest common divisor of $f$ and $g$.

```julia
lcm(f::MyElem, g::MyElem)
```

Returns $fg/gcd(f, g)$ if either $f$ or $g$ is not zero, otherwise it throws a
`DivideError()`.

```julia
gcdx(f::MyElem, g::MyElem)
```

Returns a triple `d, s, t` such that $d = gcd(f, g)$ and $d = sf + tg$, with $s$ reduced
modulo $g$ and $t$ reduced modulo $f$.

```julia
gcdinv(f::MyElem, g::MyElem)
```

Returns a tuple `d, s` such that $d = gcd(f, g)$ and $s = (f/d)^{-1} \pmod{g/d}$. Note
that $d = 1$ iff $f$ is invertible modulo $g$, in which case $s = f^{-1} \pmod{g}$.

