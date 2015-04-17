###############################################################################
#
#   pari_poly.jl : functions for Pari t_POL objects
#
###############################################################################

###############################################################################
#
#   Constructors
#
###############################################################################

PariPolyID = ObjectIdDict()

type PariPolyRing{T <: PariRing} <: PariRing
   base_ring::PariRing
   pol_0::Ptr{Int}
   S::Symbol

   function PariPolyRing(R::PariRing, s::Symbol)
      z = ccall((:pari_malloc, :libpari), Ptr{Int}, (Int,), 2*sizeof(Int))
      unsafe_store!(z, evaltyp(t_POL) | 2, 1) 
      unsafe_store!(z, evalsigne(0) | evalvarn(0), 2)
      try
         return PariPolyID[s]
      catch
         r = PariPolyID[s] = new(R, z, s)
         finalizer(r, _pari_poly_zero_clear_fn)
         return r
      end
      
   end
end

function _pari_poly_zero_clear_fn(p::PariPolyRing)
   ccall((:pari_free, :libpari), Void, (Ptr{Int},), p.pol_0)
end

type pari_poly{T <: PariRing} <: RingElem
   d::Ptr{Int}
   parent::PariPolyRing{T}

   function pari_poly(data::Ptr{Int})
      g = new(gclone(data))
      finalizer(g, _pari_poly_unclone)
      return g
   end

   function pari_poly(s::Int)
      g = new(ccall((:pari_malloc, :libpari), Ptr{Int}, (Int,), s*sizeof(Int)))
      finalizer(g, _pari_poly_clear_fn)
      return g
   end
end

function _pari_poly_clear_fn(g::pari_poly)
   ccall((:pari_free, :libpari), Void, (Ptr{Uint},), g.d)
end

_pari_poly_unclone(g::pari_poly) = gunclone(g.d)

parent{T <: PariRing}(a::pari_poly{T}) = a.parent

var(p::PariPolyRing) = p.S

###############################################################################
#
#   String I/O
#
###############################################################################

show(io::IO, x::pari_poly) = pari_print(io, x.d)

function show{T <: Ring}(io::IO, p::PariPolyRing{T})
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
   g = pari_poly{PariIntegerRing}(s)
   g.parent = PariPolyRing{PariIntegerRing}(PariZZ, var(parent(a)))
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
   c = ZZ()
   for i = 0:length - 1
      ZZ!(c, reinterpret(Ptr{Int}, unsafe_load(g, 3 + i)))
      setcoeff!(z, i, c)
   end
   return z
end

###############################################################################
#
#   Factorisation
#
###############################################################################

function factor{T <: PariRing}(pol::pari_poly{T})
   av = unsafe_load(avma, 1)
   f = ccall((:factor, :libpari), Ptr{Int}, (Ptr{Int},), pol.d)
   fac = PariFactor{PariPolyRing{T}}(f, pol.parent)
   unsafe_store!(avma, av, 1)
   return fac
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function Base.call(ord::PariPolyRing{PariIntegerRing}, n::Ptr{Int})
   pol = pari_poly{PariIntegerRing}(n)
   pol.parent = PariPolyRing{PariIntegerRing}(PariZZ)
   return pol
end


function Base.call(a::FmpzPolyRing, g::pari_poly{PariIntegerRing})
   z = fmpz_poly()
   z.parent = a
   fmpz_poly!(z, g.d)
   return z
end

