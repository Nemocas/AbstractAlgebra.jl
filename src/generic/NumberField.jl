###############################################################################
#
#   number_field.jl - Generic number fields (mainly for test code)
#
###############################################################################

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function  gen(R::ResField{PolyRingElem{Rational{BigInt}}})
   return R(gen(base_ring(R)))
end

###############################################################################
#
#   Random generation
#
###############################################################################

function RandomExtensions.make(S::ResField{Generic.Poly{Rational{BigInt}}}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1])
   else
      n = degree(S.modulus)
      make(S, make(base_ring(S), n - 1:n - 1, vs...))
   end
end

###############################################################################
#
#   number_field constructor (mainly for test code)
#
###############################################################################

function number_field(a::Poly{Rational{BigInt}}, s::Symbol, t = "\$"; cached = true)
   S = parent(a)
   R = residue_field(S, a, cached=cached)
   x = gen(S)
   return R, R(x)
end
