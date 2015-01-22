###########################################################################################
#
#   pari_frac.jl : fractions in Pari
#
###########################################################################################

###########################################################################################
#
#   Constructors
#
###########################################################################################

type PariRationalField <: PariField
end

PariQQ = PariRationalField()

type pari_rat <: RingElem
   d::Ptr{Int}

   function pari_rat(s::Int)
      g = new(ccall((:pari_malloc, :libpari), Ptr{Int}, (Int,), s*BITS_IN_WORD))
      finalizer(g, _pari_rat_clear_fn)
      return g
   end
end

_pari_rat_clear_fn(g::pari_rat) = ccall((:pari_free, :libpari), Void, (Ptr{Uint},), g.d)

parent(a::pari_rat) = PariQQ

###########################################################################################
#
#   String I/O
#
###########################################################################################

show(io :: IO, a::PariRationalField) = print(io, "Rational Field")

show(io::IO, x::pari_rat) = pari_print(io, x.d)

###########################################################################################
#
#   Conversions to from Rational{BigInt}
#
###########################################################################################

gensize(a::Rational{BigInt}) = words(a.num) + words(a.den) + 7

function pari!(x::Ptr{Int}, a::Rational{BigInt})
   if a.den == 1 # we have an integer
      return pari!(x, a.num, gensize(a.num))
   end
   unsafe_store!(x, evaltyp(t_FRAC) | 3, 1) 
   unsafe_store!(x, x + sizeof(Int)*3, 2)
   s = pari!(x + sizeof(Int)*3, a.num, gensize(a.num))
   unsafe_store!(x, x + sizeof(Int)*(s + 3), 3)
   s += pari!(x + sizeof(Int)*(3 + s), a.den, gensize(a.den))
   return s + 3
end

function pari(a::Rational{BigInt})
   s = gensize(a)
   g = pari_rat(s)
   pari!(reinterpret(Ptr{Int}, g.d), a)
   return g
end

function QQ!(z::Rational{BigInt}, g::Ptr{Int})
   if typ(g) == t_INT
      ZZ!(z.num, g)
      ZZ!(z.den, reinterpret(Ptr{Int}, unsafe_load(gen_1, 1)))
   else
      ZZ!(z.num, reinterpret(Ptr{Int}, unsafe_load(g, 2)))
      ZZ!(z.den, reinterpret(Ptr{Int}, unsafe_load(g, 3)))
   end
end

function call(::RationalField, g::pari_rat)
   r = QQ()
   QQ!(r, g.d)
   return r
end
