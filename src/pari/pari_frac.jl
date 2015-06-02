###############################################################################
#
#   pari_frac.jl : rationals and fractions in Pari
#
###############################################################################

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

PariQQ = PariRationalField()

parent(a::pari_rat) = PariQQ

###############################################################################
#
#   String I/O
#
###############################################################################

show(io :: IO, a::PariRationalField) = print(io, "Rational Field")

show(io::IO, x::pari_rat) = pari_print(io, x.d)

###############################################################################
#
#   Conversions to from Rational{BigInt}
#
###############################################################################

gensize(a::fmpq) = words(num(a)) + words(den(a)) + 7

function pari!(x::Ptr{Int}, a::fmpq)
   if den(a) == 1 # we have an integer
      return pari!(x, num(a), gensize(num(a)))
   end
   unsafe_store!(x, evaltyp(t_FRAC) | 3, 1) 
   unsafe_store!(x, x + sizeof(Int)*3, 2)
   s = pari!(x + sizeof(Int)*3, num(a), gensize(num(a)))
   unsafe_store!(x, x + sizeof(Int)*(s + 3), 3)
   s += pari!(x + sizeof(Int)*(3 + s), den(a), gensize(den(a)))
   return s + 3
end

function pari(a::fmpq)
   s = gensize(a)
   g = pari_rat(s)
   pari!(reinterpret(Ptr{Int}, g.d), a)
   return g
end

# this function should be removed once it is
# possible to pass to ZZ! in pari_int
function ZZ!(z::Ptr{fmpz}, g::Ptr{Int})
   const data_offset = 1 + div(2*sizeof(Cint), sizeof(Int))
   s = (unsafe_load(g, 1) & LGBITS) - 2
   sgn = signe(g + sizeof(Int))
   if s == 0
      ccall((:fmpz_zero, :libflint), Void, (Ptr{fmpz},), z)
   elseif s == 1
      ccall((:fmpz_set_ui, :libflint), Void, 
            (Ptr{fmpz}, Int), z, unsafe_load(g, 3))
   else
      a = ccall((:_fmpz_promote, :libflint), Ptr{BigInt}, (Ptr{fmpz},), z)
      ccall((:__gmpz_realloc2, :libgmp), Void, 
            (Ptr{BigInt}, Int), a, s*BITS_IN_WORD)
      for i in 1:s
         unsafe_store!(reinterpret(Ptr{UInt}, 
                       unsafe_load(reinterpret(Ptr{Int}, a), data_offset)), 
                       reinterpret(UInt, unsafe_load(g, i + 2)), i)
      end
      unsafe_store!(reinterpret(Ptr{Cint}, a), Cint(s), 2)
   end
   if sgn < 0
      ccall((:fmpz_neg, :libflint), Void, (Ptr{fmpz}, Ptr{fmpz}), z, z)
   end
end

function QQ!(z::fmpq, g::Ptr{Int})
   num = ccall((:fmpq_numerator_ptr, :libflint), Ptr{fmpz}, (Ptr{fmpq},), &z)
   den = ccall((:fmpq_denominator_ptr, :libflint), Ptr{fmpz}, (Ptr{fmpq},), &z)
   if typ(g) == t_INT
      ZZ!(num, g)
      ZZ!(den, reinterpret(Ptr{Int}, unsafe_load(gen_1, 1)))
   else
      ZZ!(num, reinterpret(Ptr{Int}, unsafe_load(g, 2)))
      ZZ!(den, reinterpret(Ptr{Int}, unsafe_load(g, 3)))
   end
end

function call(::RationalField, g::pari_rat)
   r = QQ()
   QQ!(r, g.d)
   return r
end
