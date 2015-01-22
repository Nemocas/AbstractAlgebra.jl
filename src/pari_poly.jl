###########################################################################################
#
#   pari_poly.jl : functions for Pari t_POL objects
#
###########################################################################################

###########################################################################################
#
#   Constructors
#
###########################################################################################

PariPolyID = ObjectIdDict()

type PariPolyRing{T <: PariRing, S} <: PariRing
   base_ring::PariRing
   pol_0::Ptr{Int}

   function PariPolyRing(R::PariRing)
      z = ccall((:pari_malloc, :libpari), Ptr{Int}, (Int,), 2*sizeof(Int))
      unsafe_store!(z, evaltyp(t_POL) | 2, 1) 
      unsafe_store!(z, evalsigne(0) | evalvarn(0), 2)
      try
         return PariPolyID[S]
      catch
         r = PariPolyID[S] = new(R, z)
         finalizer(r, _pari_poly_zero_clear_fn)
         return r
      end
      
   end
end

function _pari_poly_zero_clear_fn(p::PariPolyRing)
   ccall((:pari_free, :libpari), Void, (Ptr{Int},), p.pol_0)
end

type pari_poly{T <: PariRing, S} <: RingElem
   d::Ptr{Int}
   parent::PariPolyRing{T, S}

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

_pari_poly_clear_fn(g::pari_poly) = ccall((:pari_free, :libpari), Void, (Ptr{Uint},), g.d)

_pari_poly_unclone(g::pari_poly) = gunclone(g.d)

parent{T <: PariRing, S}(a::pari_poly{T, S}) = a.parent

###########################################################################################
#
#   String I/O
#
###########################################################################################

show(io::IO, x::pari_poly) = pari_print(io, x.d)

function show{T <: Ring, S}(io::IO, p::PariPolyRing{T, S})
   print(io, "Univariate Polynomial Ring in ")
   print(io, string(S))
   print(io, " over ")
   show(io, p.base_ring)
end

###########################################################################################
#
#   Conversions to from Poly{BigInt, S}
#
###########################################################################################

function gensize(a::fmpz_poly)
   coeffs = a.coeffs
   total = 0
   for i = 0:a.length - 1
      total += ccall((:fmpz_size, :libflint), Int, (Ptr{Int},), coeffs + i*sizeof(Int))
   end
   return total + 2*a.length + 2
end

function pari!(x::Ptr{Int}, a::fmpz_poly, s::Int)
   unsafe_store!(x, evaltyp(t_POL) | a.length + 2, 1) 
   unsafe_store!(x, evalsigne(Int(a.length != 0)) | evalvarn(0), 2)
   s = a.length + 2
   for i = 0:a.length - 1
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

function pari{S}(a::fmpz_poly{S})
   s = gensize(a)
   g = pari_poly{PariIntegerRing, S}(s)
   g.parent = PariPolyRing{PariIntegerRing, S}(PariZZ)
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
   c = BigInt()
   for i = 0:length - 1
      ZZ!(c, reinterpret(Ptr{Int}, unsafe_load(g, 3 + i)))
      setcoeff!(z, i, c)
   end
end

function Base.call{S}(a::FmpzPolyRing{BigInt, S}, g::pari_poly{PariIntegerRing, S})
   z = fmpz_poly{S}()
   z.parent = a
   fmpz_poly!(z, g.d)
   return z
end

