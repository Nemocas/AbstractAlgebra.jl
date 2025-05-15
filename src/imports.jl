using Random: Random, AbstractRNG, SamplerTrivial
using RandomExtensions: RandomExtensions, make, Make, Make2, Make3, Make4

if !isdefined(Base, :get_extension)
  using Requires: @require
end

import LinearAlgebra

import LinearAlgebra: det
import LinearAlgebra: dot
import LinearAlgebra: hessenberg
import LinearAlgebra: isdiag
import LinearAlgebra: ishermitian
import LinearAlgebra: issymmetric
import LinearAlgebra: istril
import LinearAlgebra: istriu
import LinearAlgebra: lu
import LinearAlgebra: lu!
import LinearAlgebra: norm
import LinearAlgebra: nullspace
import LinearAlgebra: rank
import LinearAlgebra: tr


# Do not import elements of `Base_import_exclude` as we have our own
import Base: !=
import Base: &
import Base: *
import Base: +
import Base: -
import Base: /
import Base: //
import Base: <
import Base: <<
import Base: <=
import Base: ==
import Base: >
import Base: >=
import Base: >>
import Base: Array
import Base: Matrix
import Base: Rational
import Base: ^
import Base: abs
import Base: axes
import Base: conj
import Base: conj!
import Base: convert
import Base: copy
import Base: deepcopy
import Base: deepcopy_internal
import Base: delete!
import Base: exponent
import Base: fill
import Base: gcd
import Base: gcdx
import Base: get
import Base: getindex
import Base: getkey
import Base: hash
import Base: hcat
import Base: intersect
import Base: invmod
import Base: isempty
import Base: isequal
import Base: isfinite
import Base: isless
import Base: isone
import Base: isqrt
import Base: isreal
import Base: iszero
import Base: iterate
import Base: lcm
import Base: length
import Base: mod
import Base: ndigits
import Base: one
import Base: parent
import Base: pop!
import Base: powermod
import Base: precision
import Base: rand
import Base: rem
import Base: reverse
import Base: reverse!
import Base: setindex!
import Base: show
import Base: similar
import Base: size
import Base: sizehint!
import Base: string
import Base: trailing_zeros
import Base: transpose
import Base: truncate
import Base: typed_hcat
import Base: typed_hvcat
import Base: vcat
import Base: xor
import Base: zero
import Base: zeros
import Base: |
import Base: ~
