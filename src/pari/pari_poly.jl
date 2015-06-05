###############################################################################
#
#   pari_poly.jl : functions for Pari t_POL objects
#
###############################################################################

export PariPolyRing, pari_poly, fmpz_poly!

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent{T <: RingElem}(a::pari_poly{T}) = a.parent

var(p::PariPolyRing) = p.S

###############################################################################
#
#   String I/O
#
###############################################################################

show(io::IO, x::pari_poly) = pari_print(io, x.d)

function show{T <: RingElem}(io::IO, p::PariPolyRing{T})
   print(io, "Univariate Polynomial Ring in ")
   print(io, string(var(p)))
   print(io, " over ")
   show(io, p.base_ring)
end

###############################################################################
#
#   Conversions to from fmpz_poly
#
###############################################################################

function gensize(a::fmpz_poly)
   coeffs = a.coeffs
   total = 0
   for i = 0:length(a) - 1
      total += ccall((:fmpz_size, :libflint), Int, 
                     (Ptr{Int},), coeffs + i*sizeof(Int))
   end
   return total + 3*length(a) + 2
end

function pari!(x::Ptr{Int}, a::fmpz_poly, s::Int)
   unsafe_store!(x, evaltyp(t_POL) | length(a) + 2, 1) 
   unsafe_store!(x, evalsigne(Int(length(a) != 0)) | evalvarn(0), 2)
   s = length(a) + 2
   for i = 0:length(a) - 1
      z = coeff(a, i)
      if z == 0
         unsafe_store!(x, unsafe_load(gen_0), i + 3)
      else
         unsafe_store!(x, x + sizeof(Int)*s, i + 3)
         s += pari!(x + sizeof(Int)*s, z, gensize(z))
      end
   end
   return s
end

function pari(a::fmpz_poly)
   s = gensize(a)
   g = pari_poly{pari_int}(s)
   g.parent = PariPolyRing{pari_int}(PariZZ, var(parent(a)))
   pari!(reinterpret(Ptr{Int}, g.d), a, s)
   return g
end

function fmpz_poly!(z::fmpz_poly, g::Ptr{Int})
   length = (unsafe_load(g, 1) & LGBITS) - 2
   fit!(z, length)
   z.length = length
   if length == 0
      return
   end
   c = fmpz()
   for i = 0:length - 1
      fmpz!(c, reinterpret(Ptr{Int}, unsafe_load(g, 3 + i)))
      setcoeff!(z, i, c)
   end
   return z
end

###############################################################################
#
#   Factorisation
#
###############################################################################

function factor{T <: RingElem}(pol::pari_poly{T})
   av = unsafe_load(avma, 1)
   f = ccall((:factor, :libpari), Ptr{Int}, (Ptr{Int},), pol.d)
   fac = PariFactor{pari_poly{T}}(f, pol.parent)
   unsafe_store!(avma, av, 1)
   return fac
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function Base.call(ord::PariPolyRing{pari_int}, n::Ptr{Int})
   pol = pari_poly{pari_int}(n)
   pol.parent = PariPolyRing{pari_int}(PariZZ, var(ord))
   return pol
end


function Base.call(a::FmpzPolyRing, g::pari_poly{pari_int})
   z = fmpz_poly()
   z.parent = a
   fmpz_poly!(z, g.d)
   return z
end

