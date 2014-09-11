###########################################################################################
#
#   Padics2.jl : Padic numbers
#
###########################################################################################    

import Rings: O, valuation, isequal

import Base: sqrt, exp, log

export O, PadicField, Padic, valuation, prime, precision, isexact, sqrt, exp, log,
       teichmuller, isequal

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

PadicBase = ObjectIdDict()

type padic_ctx
   p :: Int # can't make this a ZZ
   pinv :: Float64
   pow :: Ptr{Void}
   minpre :: Int
   maxpre :: Int
   mode :: Int
   function padic_ctx(p::ZZ)
      d = new(0, 0.0, C_NULL, 0, 0, 0)
      finalizer(d, _padic_ctx_clear_fn)
      ccall((:padic_ctx_init, :libflint), Void, (Ptr{padic_ctx}, Ptr{ZZ}, Int, Int, Cint), &d, &p, 0, 0, 1)
      return d
   end
end

function _padic_ctx_clear_fn(a :: padic_ctx)
   #ccall((:padic_ctx_clear, :libflint), Void, (Ptr{padic_ctx},), &a)
end

type Padic{S} <: Field
   u :: Int # can't make this a ZZ
   v :: Int
   N :: Int
   exact :: Bool
   function Padic()
      d = new(0, 0, 0, false)
      ccall((:padic_init2, :libflint), Void, (Ptr{Padic}, Int), &d, 0)
      finalizer(d, _Padic_clear_fn)
      return d
   end
   function Padic(a::QQ)
      n = num(a)
      sign(n) < 0 && error("Cannot create negative p-adic to infinite precision")
      z = Padic{S}()
      p = prime(Padic{S})
      d = den(a)
      if d == 1
         v = 0
      elseif d == p
         v = 1
      else
         v, r = remove(d, p)
         r != 1 && error("Denominator not a power of p when converting rational to p-adic")
      end
      ccall((:padic_set_exact_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &n, &(PadicBase[S]::padic_ctx))
      z.v -= v
      z.exact = true
      return z
   end
   function Padic(a::ZZ)
      z = Padic{S}()
      repr = bool(ccall((:padic_set_exact_fmpz, :libflint), Cint, (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), &z, &a, &(PadicBase[S]::padic_ctx)))
      !repr && error("Cannot represent negative integer exactly")
      z.exact = true
      return z
   end
   Padic(a::Int) = Padic{S}(ZZ(a))
end

function _Padic_clear_fn{S}(a :: Padic{S})
   ccall((:padic_clear, :libflint), Void, (Ptr{Padic{S}},), &a)
end

function O{S}(::Type{Padic{S}}, n::QQ)
   m = den(n)
   n <= 0 && throw(DomainError())
   if m == 1
      return O(Padic{S}, num(n))
   end
   num(n) != 1 && error("Not a power of p in p-adic O()")
   d = Padic{S}()
   if m == prime(Padic{S})
      d.N = -1
      return d
   else
     p = prime(Padic{S})
     d.N = -flog(m, p) 
     p^-d.N != m && error("Not a power of p in p-adic O()")
   end
   return d
end

function O{S}(::Type{Padic{S}}, n::ZZ)
   n <= 0 && throw(DomainError())
   d = Padic{S}()
   if n == 1
      d.N = 0
   elseif n == prime(Padic{S})
      d.N = 1
      return d
   else
     p = prime(Padic{S})
     d.N = flog(n, p) 
     p^d.N != n && error("Not a power of p in p-adic O()")
   end
   return d
end

O{S}(::Type{Padic{S}}, n::Int) = O(Padic{S}, ZZ(n))

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function prime{S}(::Type{Padic{S}})
   ctx = PadicBase[S]::padic_ctx
   z = ZZ()
   ccall((:padic_ctx_pow_ui, :libflint), Void, (Ptr{ZZ}, Int, Ptr{padic_ctx}), &z, 1, &ctx)
   return z 
end

precision{S}(a::Padic{S}) = a.N

valuation{S}(a::Padic{S}) = a.v

isexact{S}(a::Padic{S}) = a.exact

function zero{S}(::Type{Padic{S}})
   z = Padic{S}()
   ccall((:padic_zero, :libflint), Void, (Ptr{Padic},), &z)
   z.exact = true
   return z
end

function one{S}(::Type{Padic{S}})
   z = Padic{S}()
   z.N = 1
   ccall((:padic_one, :libflint), Void, (Ptr{Padic},), &z)
   z.exact = true
   return z
end

iszero{S}(a::Padic{S}) = bool(ccall((:padic_is_zero, :libflint), Cint, (Ptr{Padic},), &a)) && a.exact

isone{S}(a::Padic{S}) = bool(ccall((:padic_is_one, :libflint), Cint, (Ptr{Padic},), &a)) && a.exact

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{S}(io::IO, x::Padic{S})
   cstr = ccall((:padic_get_str, :libflint), Ptr{Uint8}, 
               (Ptr{Void}, Ptr{Padic}, Ptr{padic_ctx}), C_NULL, &x, &(PadicBase[S]::padic_ctx))

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   if !x.exact
      print(io, " + O(")
      print(io, prime(Padic{S}))
      print(io, "^$(x.N))")
   end
end

function show{S}(io::IO, ::Type{Padic{S}})
   print(io, "p-adic number field")
end

needs_parentheses{S}(x::Padic{S}) = true

is_negative{S}(x::Padic{S}) = false

show_minus_one{S}(::Type{Padic{S}}) = true

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit{S}(x::Padic{S}) = x

###########################################################################################
#
#   Unary operators
#
###########################################################################################

function -{S}(x::Padic{S})
   if iszero(x)
      return x
   end
   x.exact && error("Cannot compute infinite precision p-adic")
   z = Padic{S}()
   z.N = x.N
   ccall((:padic_neg, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &(PadicBase[S]::padic_ctx))
   return z
end

###########################################################################################
#
#   Binary operators
#
###########################################################################################

function +{S}(x::Padic{S}, y::Padic{S})
   z = Padic{S}()
   if x.exact && y.exact
      ccall((:padic_add_exact, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &y, &(PadicBase[S]::padic_ctx))
      z.exact = true
   else
      z.N = x.exact ? y.N : (y.exact ? x.N : min(x.N, y.N))
      ccall((:padic_add, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &y, &(PadicBase[S]::padic_ctx))
   end
   return z
end

function -{S}(x::Padic{S}, y::Padic{S})
   z = Padic{S}()
   if x.exact && y.exact
      repr = bool(ccall((:padic_sub_exact, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &y, &(PadicBase[S]::padic_ctx)))
      !repr && error("Unable to represent exact result of subtraction")
      z.exact = true
   else
      z.N = x.exact ? y.N : (y.exact ? x.N : min(x.N, y.N))
      ccall((:padic_sub, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &y, &(PadicBase[S]::padic_ctx))
   end
   return z
end

function *{S}(x::Padic{S}, y::Padic{S})
   z = Padic{S}()
   if x.exact && y.exact
      ccall((:padic_mul_exact, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &y, &(PadicBase[S]::padic_ctx))
      z.exact = true
   else 
      z.N = x.exact ? y.N + x.v : (y.exact ? x.N + y.v : min(x.N + y.v, y.N + x.v))
      ccall((:padic_mul, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &y, &(PadicBase[S]::padic_ctx))
   end
   return z
end

function gcd{S}(x::Padic{S}, y::Padic{S})
   if x == 0 && y == 0 
      z = Padic{S}()
      z.N = y.exact ? x.N : (x.exact ? y.N : min(x.N, y.N))
      if x.exact && y.exact
         z.exact = true
      end
   else
      z = Padic{S}(1)
   end
   return z
end

###########################################################################################
#
#   Unsafe operators
#
###########################################################################################

function mul!{S}(z::Padic{S}, x::Padic{S}, y::Padic{S})
   if x.exact && y.exact
      ccall((:padic_mul_exact, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &y, &(PadicBase[S]::padic_ctx))
      z.exact = true
   else 
      z.N = x.exact ? y.N + x.v : (y.exact ? x.N + y.v : min(x.N + y.v, y.N + x.v))
      ccall((:padic_mul, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &y, &(PadicBase[S]::padic_ctx))
   end
end

function addeq!{S}(x::Padic{S}, y::Padic{S})
   if x.exact && y.exact
      ccall((:padic_add_exact, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &x, &x, &y, &(PadicBase[S]::padic_ctx))
      x.exact = true
   else
      x.N = x.exact ? y.N : (y.exact ? x.N : min(x.N, y.N))
      ccall((:padic_add, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &x, &x, &y, &(PadicBase[S]::padic_ctx))
   end
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function +{S}(x::Padic{S}, y::QQ)
   n = num(y)
   if sign(n) < 0
      return x - (-y)
   end
   z = Padic{S}()
   if x.exact
      p = prime(Padic{S})
      d = den(y)
      if d == 1
         v = 0
      elseif d == p
         v = 1
      else
         v, r = remove(d, p)
         r != 1 && error("Denominator not a power of p in expression involving p-adic")
      end
      ccall((:padic_set_exact_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &n, &(PadicBase[S]::padic_ctx))
      z.v -= v
      ccall((:padic_add_exact, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &x, &(PadicBase[S]::padic_ctx))
      z.exact = true
   else
      z.N = x.N
      ccall((:padic_set_fmpq, :libflint), Void, 
                (Ptr{Padic}, Ptr{Fraction}, Ptr{padic_ctx}), 
               &z, &y, &(PadicBase[S]::padic_ctx))
      ccall((:padic_add, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &x, &(PadicBase[S]::padic_ctx))
   end
   return z
end

+{S}(x::QQ, y::Padic{S}) = y + x

function +{S}(x::Padic{S}, y::ZZ)
   if sign(y) < 0
      return x - (-y)
   end
   z = Padic{S}()
   ccall((:padic_set_exact_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &y, &(PadicBase[S]::padic_ctx))
   if x.exact
      ccall((:padic_add_exact, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &x, &(PadicBase[S]::padic_ctx))
      z.exact = true
   else
      z.N = x.N
      ccall((:padic_add, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &x, &(PadicBase[S]::padic_ctx))
   end
   return z
end

+{S}(x::Padic{S}, y::Int) = x + ZZ(y)

+{S}(x::ZZ, y::Padic{S}) = y + x

+{S}(x::Int, y::Padic{S}) = y + ZZ(x)

function -{S}(x::Padic{S}, y::QQ)
   n = num(y)
   if sign(n) < 0
      return x - (-y)
   end
   z = Padic{S}()
   if x.exact
      p = prime(Padic{S})
      d = den(y)
      if d == 1
         v = 0
      elseif d == p
         v = 1
      else
         v, r = remove(d, p)
         r != 1 && error("Denominator not a power of p in expression involving p-adic")
      end
      ccall((:padic_set_exact_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &n, &(PadicBase[S]::padic_ctx))
      z.v -= v
      ccall((:padic_sub_exact, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &z, &(PadicBase[S]::padic_ctx))
      z.exact = true
   else
      z.N = x.N
      ccall((:padic_set_fmpq, :libflint), Void, 
                (Ptr{Padic}, Ptr{QQ}, Ptr{padic_ctx}), 
               &z, &y, &(PadicBase[S]::padic_ctx))
      ccall((:padic_sub, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &z, &(PadicBase[S]::padic_ctx))
   end
   return z
end

function -{S}(x::Padic{S}, y::ZZ)
   if sign(y) < 0
      return x + (-y)
   end
   z = Padic{S}()
   ccall((:padic_set_exact_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &y, &(PadicBase[S]::padic_ctx))
   if x.exact
      repr = bool(ccall((:padic_sub_exact, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &z, &(PadicBase[S]::padic_ctx)))
      !repr && error("Unable to represent negative value exactly")
      z.exact = true
   else
      z.N = x.N
      ccall((:padic_sub, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &x, &z, &(PadicBase[S]::padic_ctx))
   end
   return z
end

function -{S}(x::QQ, y::Padic{S})
   n = num(x)
   if sign(n) < 0
      return -y + (-x)
   end
   z = Padic{S}()
   if y.exact
      p = prime(Padic{S})
      d = den(x)
      if d == 1
         v = 0
      elseif d == p
         v = 1
      else
         v, r = remove(d, p)
         r != 1 && error("Denominator not a power of p in expression involving p-adic")
      end
      ccall((:padic_set_exact_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &n, &(PadicBase[S]::padic_ctx))
      z.v -= v
      repr = bool(ccall((:padic_sub_exact, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &y, &(PadicBase[S]::padic_ctx)))
      !repr && error("Unable to represent negative value exactly")
      z.exact = true
   else
      z.N = y.N
      ccall((:padic_set_fmpq, :libflint), Void, 
                (Ptr{Padic}, Ptr{Fraction}, Ptr{padic_ctx}), 
               &z, &x, &(PadicBase[S]::padic_ctx))
      ccall((:padic_sub, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &y, &(PadicBase[S]::padic_ctx))
   end
   return z
end

function -{S}(x::ZZ, y::Padic{S})
   if sign(x) < 0
      return -y + (-x)
   end
   z = Padic{S}()
   ccall((:padic_set_exact_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &x, &(PadicBase[S]::padic_ctx))
   if y.exact
      repr = bool(ccall((:padic_sub_exact, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &y, &(PadicBase[S]::padic_ctx)))
      !repr && error("Unable to represent negative value exactly")
      z.exact = true
   else
      z.N = y.N
      ccall((:padic_sub, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &y, &(PadicBase[S]::padic_ctx))
   end
   return z
end

-{S}(x::Padic{S}, y::Int) = x - ZZ(y)

-{S}(x::Int, y::Padic{S}) = ZZ(x) - y

function *{S}(x::Padic{S}, y::QQ)
   n = num(y)
   if sign(n) < 0
      return -(x*(-y))
   end
   if y == 0
      return zero(Padic{S})
   end
   z = Padic{S}()
   d = den(y)
   if x.exact
      z.exact = true
      res = bool(ccall((:padic_div_exact_fmpz, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &x, &d, &(PadicBase[S]::padic_ctx)))
      !res && error("Unable to multiply by rational to infinite precision")
      c = Padic{S}()
      ccall((:padic_set_exact_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &c, &n, &(PadicBase[S]::padic_ctx))
      ccall((:padic_mul_exact, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &c, &(PadicBase[S]::padic_ctx))
   else
      p = prime(Padic{S})
      v1, r = remove(n, p)
      v2, r = remove(d, p)
      z.N = x.N - x.v + v1 - v2
      ccall((:padic_set_fmpq, :libflint), Void, 
                (Ptr{Padic}, Ptr{Fraction}, Ptr{padic_ctx}), 
               &z, &y, &(PadicBase[S]::padic_ctx))
      z.N = x.N + z.v
      ccall((:padic_mul, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &x, &(PadicBase[S]::padic_ctx))
   end
   return z
end

*{S}(x::QQ, y::Padic{S}) = y*x

function *{S}(x::Padic{S}, y::ZZ)
   if sign(y) < 0
      return -(x*(-y))
   end
   if y == 0
      return zero(Padic{S})
   end
   z = Padic{S}()
   ccall((:padic_set_exact_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &y, &(PadicBase[S]::padic_ctx))
   if x.exact
      z.exact = true
      ccall((:padic_mul_exact, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &x, &(PadicBase[S]::padic_ctx))
   else
      z.N = x.N + z.v
      ccall((:padic_mul, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &z, &x, &(PadicBase[S]::padic_ctx))
   end
   return z
end

*{S}(x::ZZ, y::Padic{S}) = y*x

*{S}(x::Padic{S}, y::Int) = x*ZZ(y)

*{S}(x::Int, y::Padic{S}) = y*ZZ(x)

###########################################################################################
#
#   Comparison
#
###########################################################################################

function =={S}(a::Padic{S}, b::Padic{S})
   if a.exact && b.exact
      return bool(ccall((:padic_equal, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}), &a, &b))
   else
      z = Padic{S}()
      z.N = a.exact ? b.N : (b.exact ? a.N : min(a.N, b.N))
      ccall((:padic_sub, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), 
               &z, &a, &b, &(PadicBase[S]::padic_ctx))
      return bool(ccall((:padic_is_zero, :libflint), Cint, 
                (Ptr{Padic},), &z))
   end
end

function isequal{S}(a::Padic{S}, b::Padic{S})
   if a.exact != b.exact
      return false
   end
   if a.exact == true
      return a == b
   else
      return a.N == b.N && a == b
   end
end

###########################################################################################
#
#   Ad hoc comparison
#
###########################################################################################

function =={S}(a::Padic{S}, b::QQ)
   if a.exact
      if sign(num(b)) < 0
         return false
      end
      z = Padic{S}(b)
   else
      z = Padic{S}()
      z.N = a.N
      ccall((:padic_set_fmpq, :libflint), Void, 
                (Ptr{Padic}, Ptr{Fraction}, Ptr{padic_ctx}), 
               &z, &b, &(PadicBase[S]::padic_ctx))
   end
   return bool(ccall((:padic_equal, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}), &a, &z))
end

function =={S}(a::Padic{S}, b::ZZ)
   if a.exact
      if sign(b) < 0
         return false
      end
      z = Padic{S}(b)
   else
      z = Padic{S}()
      z.N = a.N
      ccall((:padic_set_fmpz, :libflint), Void, 
                (Ptr{Padic}, Ptr{ZZ}, Ptr{padic_ctx}), 
               &z, &b, &(PadicBase[S]::padic_ctx))
   end
   return bool(ccall((:padic_equal, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}), &a, &z))
end

=={S}(a::Padic{S}, b::Int) = a == ZZ(b)

=={S}(a::ZZ, b::Padic{S}) = b == a

=={S}(a::Int, b::Padic{S}) = b == ZZ(a)

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{S}(a::Padic{S}, n::Int)
   z = Padic{S}()
   if a.exact
      repr = bool(ccall((:padic_pow_exact_si, :libflint), Cint, 
                (Ptr{Padic}, Ptr{Padic}, Int, Ptr{padic_ctx}), 
               &z, &a, n, &(PadicBase[S]::padic_ctx)))
      !repr && error("Unable to invert p-adic to infinite precision")
      z.exact = true
   else
      z.N = a.N + (n - 1)*a.v
      ccall((:padic_pow_si, :libflint), Void, 
                (Ptr{Padic}, Ptr{Padic}, Int, Ptr{padic_ctx}), 
               &z, &a, n, &(PadicBase[S]::padic_ctx))
   end
   return z
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S}(a::Padic{S}, b::Padic{S})
   b == 0 && throw(DivideError())
   z = Padic{S}()
   if a.exact
      if b.exact
         z.exact = true
         repr = bool(ccall((:padic_div_exact, :libflint), Cint, (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}), &z, &a, &b))
         !repr && error("Unable to compute quotient of p-adics to infinite precision")
         return z
      end
      z.N = b.N - 2*b.v + a.v
   elseif b.exact
      z.N = a.N - b.v
   else
      z.N = min(a.N - b.v, b.N - 2*b.v + a.v)
   end
   ccall((:padic_div, :libflint), Cint, (Ptr{Padic}, Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), &z, &a, &b, &(PadicBase[S]::padic_ctx))
   return z
end

/{S}(a::Padic{S}, b::Padic{S}) = divexact(a, b)

divexact{S}(a::Padic{S}, b::Int) = a*(ZZ(1)/ZZ(b))

divexact{S}(a::Padic{S}, b::ZZ) = a*(1/b)

divexact{S}(a::Padic{S}, b::QQ) = a*inv(b)

/{S}(a::Padic{S}, b::Int) = divexact(a, ZZ(b))

/{S}(a::Padic{S}, b::ZZ) = divexact(a, b)

/{S}(a::Padic{S}, b::QQ) = divexact(a, b)

divexact{S}(a::Int, b::Padic{S}) = ZZ(a)*inv(b)

divexact{S}(a::ZZ, b::Padic{S}) = inv((ZZ(1)/a)*b)

divexact{S}(a::QQ, b::Padic{S}) = inv(inv(a)*b)

/{S}(a::Int, b::Padic{S}) = divexact(ZZ(a), b)

/{S}(a::ZZ, b::Padic{S}) = divexact(a, b)

/{S}(a::QQ, b::Padic{S}) = divexact(a, b)

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{S}(a::Padic{S})
   a == 0 && throw(DivideError())
   z = Padic{S}()
   if a.exact
      z.exact = true
      repr = bool(ccall((:padic_inv_exact, :libflint), Cint, (Ptr{Padic}, Ptr{Padic}), &z, &a))      
      !repr && error("Unable to invert infinite precision p-adic")
   else
      z.N = a.N - 2*a.v
      ccall((:padic_inv, :libflint), Cint, (Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), &z, &a, &(PadicBase[S]::padic_ctx))
   end
   return z
end

###########################################################################################
#
#   Square root
#
###########################################################################################

function sqrt{S}(a::Padic{S})
   (a.v % 2) != 0 && error("Unable to take padic square root")
   z = Padic{S}()
   if a.exact
      res = bool(ccall((:padic_sqrt_exact, :libflint), Cint, (Ptr{Padic}, Ptr{Padic}), &z, &a))      
      !res && error("Unable to take square root of p-adic to infinite precision")
      z.exact = true
      return z
   end
   z.N = a.N - div(a.v, 2)
   res = bool(ccall((:padic_sqrt, :libflint), Cint, (Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), &z, &a, &(PadicBase[S]::padic_ctx)))      
   !res && error("Square root of p-adic does not exist")
   return z
end

###########################################################################################
#
#   Special functions
#
###########################################################################################

function exp{S}(a::Padic{S}) 
   a != 0 && a.v <= 0 && throw(DomainError())
   if a.exact
      if iszero(a)
         return one(Padic{S})
      end
      error("Unable to compute exponential of infinite precision value")
   end
   z = Padic{S}()
   z.N = a.N
   res = bool(ccall((:padic_exp, :libflint), Cint, (Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), &z, &a, &(PadicBase[S]::padic_ctx)))
   !res && error("Unable to compute exponential")
   return z
end

function log{S}(a::Padic{S}) 
   (a.v > 0 || a.v < 0 || a == 0) && throw(DomainError())
   if a.exact
      if isone(a)
         return zero(Padic{S})
      end
      error("Unable to compute logarithm of infinite precision value")
   end
   z = Padic{S}()
   z.N = a.N
   res = bool(ccall((:padic_log, :libflint), Cint, (Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), &z, &a, &(PadicBase[S]::padic_ctx)))
   !res && error("Unable to compute logarithm")
   return z
end

function teichmuller{S}(a::Padic{S})
   a.v < 0 && throw(DomainError())
   if a.exact
      if isone(a) || iszero(a)
         return a
      end
      error("Unable to compute Teichmuller unit for infinite precision value")
   end
   z = Padic{S}()
   z.N = a.N
   ccall((:padic_teichmuller, :libflint), Void, (Ptr{Padic}, Ptr{Padic}, Ptr{padic_ctx}), &z, &a, &(PadicBase[S]::padic_ctx))
   return z
end
  
###########################################################################################
#
#   Conversions and promotions
#
###########################################################################################

convert{S}(::Type{Padic{S}}, x::Int) = Padic{S}(x)

convert{S}(::Type{Padic{S}}, x::ZZ) = Padic{S}(x)

convert{S}(::Type{Padic{S}}, x::QQ) = Padic{S}(x)

promote_rule{S}(::Type{Padic{S}}, ::Type{Int}) = Padic{S}

promote_rule{S}(::Type{Padic{S}}, ::Type{ZZ}) = Padic{S}

promote_rule{S}(::Type{Padic{S}}, ::Type{QQ}) = Padic{S}

###########################################################################################
#
#   PadicField constructor
#
###########################################################################################

function PadicField(p::ZZ)
   !isprime(p) && error("Prime base required in PadiCNumberField")
   S = gensym("padic")
   R = Padic{S}
   PadicBase[S] = padic_ctx(p)
   return R
end

PadicField(p::Int) = PadicField(ZZ(p))
