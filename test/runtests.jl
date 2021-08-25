using AbstractAlgebra

using Test

R, x = PolynomialRing(QQ, "x")

a = x//(x + 1)

b = QQ(1//2)

S, a = RationalFunctionField(QQ, "a")

y = S(a)

y*b

