```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = :(using AbstractAlgebra)
```

# Euclidean Ring Interface

If a ring provides a meaningful Euclidean structure such that a useful Euclidean
remainder can be computed practically, various additional functionality is provided
by AbstractAlgebra.jl for those rings. This functionality depends on the following
functions existing. An implementation must provide `divrem`, and the remaining
are optional as generic fallbacks exist.

```@docs
divrem
mod(f::T, g::T) where T <: RingElem
Base.div(f::T, g::T) where T <: RingElem
mulmod(f::T, g::T, m::T) where T <: RingElem
powermod(f::T, e::Int, m::T) where T <: RingElem
invmod(f::T, m::T) where T <: RingElem
divides(f::T, g::T) where T <: RingElem
remove(f::T, p::T) where T <: RingElem
valuation(f::T, p::T) where T <: RingElem
gcd(f::T, g::T) where T <: RingElem
gcd(f::T, g::T, hs::T...) where T <: RingElem
gcd(fs::AbstractArray{<:T}) where T <: RingElem
lcm(f::T, g::T) where T <: RingElem
lcm(f::T, g::T, hs::T...) where T <: RingElem
lcm(fs::AbstractArray{<:T}) where T <: RingElem
gcdx(f::T, g::T) where T <: RingElem
gcdinv(f::T, g::T) where T <: RingElem
crt(r1::T, m1::T, r1::T, m2::T; check::Bool=true) where T <: RingElement
crt(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
crt_with_lcm(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
crt_with_lcm(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
```
