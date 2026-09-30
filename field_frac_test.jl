

# [Oliver setup]
#
# activate ~/oscar-puiseux
#



R, (t,) = puiseux_polynomial_ring(QQ, ["t"])

F = fraction_field(R)
a = F(t^(1//2)) // F(t + 1)
b = (t - 1) // (t^(1//2) - 1)  # this lies in R but it's a conversion to fraction field using //
typeof(b)
c = 1//(1+t)
d = t+t^2
typeof(d)


# No difference for multivariate case
S, (x, y) = puiseux_polynomial_ring(QQ, ["x", "y"])
G = fraction_field(S)

e = 1 // (x + y)
f = (x^2 + y^2) // (x - y)
e*f
e+f


# Matrices over fraction fields

R, (s, t) = puiseux_polynomial_ring(QQ, ["s", "t"])
F = fraction_field(R)

A = matrix(F, [[s^(1//2), t],[1, s + t]])
B = matrix(F, [[s^(1//2), t, 1], [1, s + t, 0], [s^(1//2) + 1, s + 2*t, 1]])
b = matrix(F, 2, 1, [1, t^(1//3)])
Fx, x = polynomial_ring(F, "x")

det(A)
rank(A)

det(B)
rank(B)

K = kernel(B; side = :right)
is_zero(B * K)

charpoly(Fx, A)
charpoly(Fx, B)

# need divides(::PuiseuxMPolyRingElem, ::PuiseuxMPolyRingElem)
Ainv = inv(A)
A * Ainv
A * Ainv == identity_matrix(F, 2)

y = solve(A, b; side = :right)
A * y
A * y == b


############################################
#
# Stuff in Laurent that is not in Puiseux
#
############################################

`From claude:`
I compared the methods defined for Laurent polynomials in this repo with those for Puiseux, then tried each missing one on a Puiseux element:

    `divides`
    > Needed for `inv` and `solve`.

    `leading_coefficient`, `leading_term`, `leading_monomial`, `leading_exponent_vector`
    > Easy to delegate to `poly(f)`; the exponent vector needs dividing by `scale(f)`. `is_unit` already uses `leading_coefficient(poly(f))` internally.

    `terms`, `iterate`, `coeff(f, i)`, `coeff(f, m)`
    > `coefficients`, `monomials` and `exponent_vectors` exist, but there's no iteration over terms or indexed coefficient access.

    `gen(R, i)`, `var_index`
    > Only `gens(R)` exists.

    `evaluate`
    > Needs a choice: evaluating `t^(1//2)` at a point needs roots in the target ring. Evaluating at Puiseux elements, or at `x^N`-type values, is well-defined.

    `derivative`
    > Well-defined: d/dt of t^(p/q) is (p/q)·t^(p/q − 1).

    `map_coefficients`, `change_base_ring`
    > Straightforward delegation, but the constructor currently requires a `Field`.

    `factor`, `factor_squarefree`
    > **Not meaningful in the usual sense.** The ring isn't a UFD because `t - 1` keeps splitting forever. At best you'd get factorisation at a fixed scale.

    `inflate`
    > The code has comments saying to switch `rescale`/`divexact` to `inflate` once available. It's in AbstractAlgebra now, so that's a free cleanup.

    `rand`, `push_term!`/`MPolyBuildCtx`
    > Utilities for building elements and for tests.

    `add!`, `sub!`, `neg!`, `zero!`
    > In-place operations. Performance only, but matrix code calls them a lot.
