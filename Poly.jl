###########################################################################################
#
#   Poly.jl : Polynomials over rings
#
###########################################################################################    

export Poly, PolynomialRing, coeff, zero, one, gen, isgen, normalise, chebyshev_t,
       chebyshev_u, theta_qexp, eta_qexp, swinnerton_dyer, cos_minpoly, cyclotomic,
       pseudorem, pseudodivrem, primpart, content, divexact, evaluate, compose, deriv,
       resultant, lead, discriminant, bezout, truncate, mullow, divrem, mulmod, powmod,
       invmod, canonical_unit, integral, lcm, reverse, shift_left, shift_right

import Base: convert, zero, show

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

type fmpz_poly <: Ring
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   function fmpz_poly(coeffs::Ptr{Void}, alloc::Int, length::Int)
      d = new(coeffs, alloc, length)
      finalizer(d, _fmpz_poly_clear_fn)
      return d
   end
end

type fmpz_mod_poly <: Ring
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   p :: Int # can't make this a ZZ
   function fmpz_mod_poly(coeffs::Ptr{Void}, alloc::Int, length::Int, p::Int)
      d = new(coeffs, alloc, length, p)
      finalizer(d, _fmpz_mod_poly_clear_fn)
      return d
   end  
end

type fmpq_poly <: Ring
   coeffs :: Ptr{Void}
   den :: Int # can't make this a ZZ
   alloc :: Int
   length :: Int
   function fmpq_poly(coeffs::Ptr{Void}, den::Int, alloc::Int, length::Int)
      d = new(coeffs, den, alloc, length)
      finalizer(d, _fmpq_poly_clear_fn)
      return d
   end  
end

type PolyStruct{T <: Ring}
   coeffs :: Array{T, 1}
   length :: Int
end

type Poly{T <: Ring, S} <: Ring
   data :: Union(fmpz_poly, fmpz_mod_poly, fmpq_poly, PolyStruct{T})
   
   Poly(a :: PolyStruct{T}) = new(a)   

   Poly(a::fmpz_mod_poly) = new(a)
   
   Poly(a :: fmpz_poly) = new(a)
   
   Poly(a :: fmpq_poly) = new(a)
   
   Poly() = Poly(Poly{T, S}, Array(T, 0))
   Poly(a::Integer) = a == 0 ? Poly(Poly{T, S}, Array(T, 0)) : Poly(Poly{T, S}, [T(a)])
   Poly(a::T) = Poly(Poly{T, S}, [a])
   Poly(a::Poly{T, S}) = a
   Poly{R <: Ring}(a::R) = convert(Poly{T, S}, a)
end

function Poly{T, S}(::Type{Poly{T, S}}, a :: Array{T, 1})
   len = length(a)
   d = PolyStruct(a, len)
   z = Poly{T, S}(d)
   z.data.length = normalise(z, len)
   return z
end

function Poly{S}(::Type{Poly{ZZ, S}}, a :: Array{ZZ, 1})
   z = fmpz_poly(C_NULL, 0, 0)
   ccall((:fmpz_poly_init2, :libflint), Void, (Ptr{fmpz_poly}, Int), &z, length(a))
   for i = 1:length(a)
      ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Void, (Ptr{fmpz_poly}, Int, Ptr{ZZ}),
         &z, i - 1, &a[i])
   end
   return Poly{ZZ, S}(z)
end   

function Poly{M, S}(::Type{Poly{Residue{ZZ, M}, S}}, a :: Array{Residue{ZZ, M}, 1})
   z = fmpz_mod_poly(C_NULL, 0, 0, 0)
   m = modulus(Residue{ZZ, M})
   ccall((:fmpz_mod_poly_init2, :libflint), Void, (Ptr{fmpz_mod_poly}, Ptr{ZZ}, Int), &z, &m, length(a))
   for i = 1:length(a)
      ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Void, (Ptr{fmpz_mod_poly}, Int, Ptr{ZZ}),
            &z, i - 1, &(a[i].data))
   end
   return Poly{Residue{ZZ, M}, S}(z)
end

function _fmpz_poly_clear_fn(a :: fmpz_poly)
   ccall((:fmpz_poly_clear, :libflint), Void, (Ptr{fmpz_poly},), &a)
end
   
function _fmpz_mod_poly_clear_fn(a :: fmpz_mod_poly)
   ccall((:fmpz_mod_poly_clear, :libflint), Void, (Ptr{fmpz_mod_poly},), &a)
end
   
function _fmpq_poly_clear_fn(a :: fmpq_poly)
   ccall((:fmpq_poly_clear, :libflint), Void, (Ptr{fmpq_poly},), &a)
end
   
###########################################################################################
#
#   Basic manipulation
#
###########################################################################################    
   
function normalise{T <: Ring, S}(a::Poly{T, S}, len::Int)
   while len > 0 && a.data.coeffs[len] == 0 # cannot use coeff(a, len - 1) here
      len -= 1
   end

   return len
end

function coeff{S}(x::Poly{ZZ, S}, n::Int)
   z = ZZ()
   ccall((:fmpz_poly_get_coeff_fmpz, :libflint), Void, (Ptr{ZZ}, Ptr{fmpz_poly}, Int), &z, &(x.data), n)
   return z
end

function coeff{S, M}(x::Poly{Residue{ZZ, M}, S}, n::Int)
   z = ZZ()
   ccall((:fmpz_mod_poly_get_coeff_fmpz, :libflint), Void, (Ptr{ZZ}, Ptr{fmpz_mod_poly}, Int), &z, &(x.data), n)
   return Residue{ZZ, M}(z)
end

coeff{T <: Ring, S}(a::Poly{T, S}, n::Int) = n >= a.data.length ? 0 : a.data.coeffs[n + 1]

lead{T <: Ring, S}(a::Poly{T, S}) = a.data.length == 0 ? zero(T) : coeff(a, a.data.length - 1)

isgen{S}(x::Poly{ZZ, S}) = bool(ccall((:fmpz_poly_is_x, :libflint), Int, (Ptr{fmpz_poly},), &(x.data)))

isgen{S, M}(x::Poly{Residue{ZZ, M}, S}) = bool(ccall((:fmpz_mod_poly_is_x, :libflint), Int, (Ptr{fmpz_mod_poly},), &(x.data)))

isgen{T <: Ring, S}(a::Poly{T, S}) = a.data.length == 2 && a.data.coeffs[1] == 0 && a.data.coeffs[2] == 1

zero{T <: Ring, S}(::Type{Poly{T, S}}) = Poly{T, S}(0)

one{T <: Ring, S}(::Type{Poly{T, S}}) = Poly{T, S}(1)

gen{T <: Ring, S}(::Type{Poly{T, S}}) = Poly(Poly{T, S}, [T(0), T(1)])

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{S}(io::IO, x::Poly{ZZ, S})
   if x.data.length == 0
      print(io, "0")
   else
      cstr = ccall((:fmpz_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{fmpz_poly}, Ptr{Uint8}), &(x.data), bytestring(string(S)))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
end

function show{S, M}(io::IO, x::Poly{Residue{ZZ, M}, S})
   if x.data.length == 0
      print(io, "0")
   else
      cstr = ccall((:fmpz_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{fmpz_mod_poly}, Ptr{Uint8}), &(x.data), bytestring(string(S)))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
end

function show{T <: Ring, S}(io::IO, x::Poly{T, S})
   len = x.data.length

   if len == 0
      print(io, zero(T))
   else
      for i = 1:len - 1
         c = x.data.coeffs[len - i + 1]
         bracket = needs_parentheses(c)
         if c != 0
            if i != 1 && !is_negative(c)
               print(io, "+")
            end
            if c != 1 && (c != -1 || show_minus_one(typeof(c)))
               if bracket
                  print(io, "(")
               end
               show(io, c)
               if bracket
                  print(io, ")")
               end
               print(io, "*")
            end
            if c == -1
               print(io, "-")
            end
            print(io, string(S))
            if len - i != 1
               print(io, "^")
               print(io, len - i)
            end
         end
      end
      c = x.data.coeffs[1]
      bracket = needs_parentheses(c)
      if c != 0
         if len != 1 && !is_negative(c)
            print(io, "+")
         end
         if bracket
            print(io, "(")
         end
         show(io, c)
         if bracket
            print(io, ")")
         end
      end
   end
end

function show{T <: Ring, S}(io::IO, ::Type{Poly{T, S}})
   print(io, "Univariate polynomial ring in ", string(S), " over ")
   show(io, T)
end

needs_parentheses{T <: Ring, S}(x::Poly{T, S}) = x.data.length > 1

is_negative{T <: Ring, S}(x::Poly{T, S}) = x.data.length <= 1 && is_negative(coeff(x, 0))

show_minus_one{T <: Ring, S}(::Type{Poly{T, S}}) = show_minus_one(T)

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit{T <: Ring, S}(x::Poly{T, S}) = canonical_unit(lead(x))

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -{S}(x::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_neg, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data))
   return z
end

function -{S, M}(x::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_neg, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(z.data), &(x.data))
   return z
end

function -{T <: Ring, S}(a::Poly{T, S})
   len = a.data.length
   d = Array(T, len)
   for i = 1:len
      d[i] = -a.data.coeffs[i]
   end
   z = Poly(Poly{T, S}, d)
   z.data.length = len
   return z
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

function +{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_add, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function -{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_sub, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function *{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_mul, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function +{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_add, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function -{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_sub, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function *{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_mul, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function +{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   lena = a.data.length
   lenb = b.data.length
   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = a.data.coeffs[i] + b.data.coeffs[i]
      i += 1
   end

   while i <= lena
      d[i] = a.data.coeffs[i]
      i += 1
   end

   while i <= lenb
      d[i] = b.data.coeffs[i]
      i += 1
   end

   z = Poly(Poly{T, S}, d)

   z.data.length = normalise(z, i - 1)

   return z
end
  
function -{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   lena = a.data.length
   lenb = b.data.length
   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = a.data.coeffs[i] - b.data.coeffs[i]
      i += 1
   end

   while i <= lena
      d[i] = a.data.coeffs[i]
      i += 1
   end

   while i <= lenb
      d[i] = -b.data.coeffs[i]
      i += 1
   end

   z = Poly(Poly{T, S}, d)

   z.data.length = normalise(z, i - 1)

   return z
end

function *{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   lena = a.data.length
   lenb = b.data.length

   if lena == 0 || lenb == 0
      return Poly{T, S}()
   end

   t = T()

   lenz = lena + lenb - 1
   d = Array(T, lenz)

   for i = 1:lena
      d[i] = a.data.coeffs[i]*b.data.coeffs[1]
   end

   for i = 2:lenb
      d[lena + i - 1] = a.data.coeffs[lena]*b.data.coeffs[i]
   end
   
   for i = 1:lena - 1
      for j = 2:lenb
         mul!(t, a.data.coeffs[i], b.data.coeffs[j])
         addeq!(d[i + j - 1], t)
      end
   end
   
   z = Poly(Poly{T, S}, d)
        
   z.data.length = normalise(z, lenz)

   return z
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function fit!{T <: Ring, S}(c::Poly{T, S}, n::Int)
   if c.data.length < n
      t = c.data.coeffs
      c.data.coeffs = Array(T, n)
      for i = 1:c.data.length
         c.data.coeffs[i] = t[i]
      end
      for i = c.data.length + 1:n
         c.data.coeffs[i] = zero(T)
      end
   end
end

function setcoeff!{T <: Ring, S}(c::Poly{T, S}, n::Int, a::T)
   if a != 0 || n + 1 <= c.data.length
      fit!(c, n + 1)
      c.data.coeffs[n + 1] = a
      c.data.length = max(c.data.length, n + 1)
      # don't normalise
   end
end

function setcoeff!{S}(z::Poly{ZZ, S}, n::Int, x::ZZ)
   ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Void, 
                (Ptr{fmpz_poly}, Int, Ptr{ZZ}), 
               &(z.data), n, &x)
end

function setcoeff!{S, M}(z::Poly{Residue{ZZ, M}, S}, n::Int, x::Residue{ZZ, M})
   ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Int, Ptr{ZZ}), 
               &(z.data), n, &(x.data))
end

function mul!{T <: Ring, S}(c::Poly{T, S}, a::Poly{T, S}, b::Poly{T, S})
   lena = a.data.length
   lenb = b.data.length

   if lena == 0 || lenb == 0
      c.data.length = 0
   else
      t = T()

      lenc = lena + lenb - 1
      fit!(c, lenc)

      for i = 1:lena
         mul!(c.data.coeffs[i], a.data.coeffs[i], b.data.coeffs[1])
      end

      for i = 2:lenb
         mul!(c.data.coeffs[lena + i - 1], a.data.coeffs[lena], b.data.coeffs[i])
      end

      for i = 1:lena - 1
         for j = 2:lenb
            mul!(t, a.data.coeffs[i], b.data.coeffs[j])
            addeq!(c.data.coeffs[i + j - 1], t)
         end
      end
        
      c.data.length = normalise(c, lenc)
   end
end

function addeq!{T <: Ring, S}(c::Poly{T, S}, a::Poly{T, S})
   lenc = c.data.length
   lena = a.data.length
   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      addeq!(c.data.coeffs[i], a.data.coeffs[i])
   end
   c.data.length = normalise(c, len)
end

function mul!{S}(z::Poly{ZZ, S}, x::Poly{ZZ, S}, y::Poly{ZZ, S})
   ccall((:fmpz_poly_mul, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data), &(y.data))
end

function mul!{S, M}(z::Poly{Residue{ZZ, M}, S}, x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   ccall((:fmpz_mod_poly_mul, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(z.data), &(x.data), &(y.data))
end

function addeq!{S}(z::Poly{ZZ, S}, x::Poly{ZZ, S},)
   ccall((:fmpz_poly_add, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(z.data), &(x.data))
end

function addeq!{S, M}(z::Poly{Residue{ZZ, M}, S}, x::Poly{Residue{ZZ, M}, S},)
   ccall((:fmpz_mod_poly_add, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(z.data), &(z.data), &(x.data))
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{S}(x::Int, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_scalar_mul_si, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), 
               &(z.data), &(y.data), x)
   return z
end

function *{S}(x::ZZ, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{ZZ}), 
               &(z.data), &(y.data), &x)
   return z
end

function *{S, M}(x::ZZ, y::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{ZZ}), 
               &(z.data), &(y.data), &x)
   return z
end

*{S, M}(x::Int, y::Poly{Residue{ZZ, M}, S}) = ZZ(x)*y

*{S, M}(x::Residue{ZZ, M}, y::Poly{Residue{ZZ, M}, S}) = x.data*y

function *{T <: Ring, S}(a::Int, b::Poly{T, S})
   len = b.data.length
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = Poly(Poly{T, S}, d)
   z.data.length = normalise(z, len)
   return z
end

function *{T <: Ring, S}(a::ZZ, b::Poly{T, S})
   len = b.data.length
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = Poly(Poly{T, S}, d)
   z.data.length = normalise(z, len)
   return z
end

function +{S}(x::Poly{ZZ, S}, y::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_add_si, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

function +{S}(x::Poly{ZZ, S}, y::ZZ)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_add_fmpz, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{ZZ}), 
               &(z.data), &(x.data), &y)
   return z
end

function +{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Int)
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_add_si, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

function +{S, M}(x::Poly{Residue{ZZ, M}, S}, y::ZZ)
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_add_fmpz, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{ZZ}), 
               &(z.data), &(x.data), &y)
   return z
end

+{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Residue{ZZ, M}) = x + y.data

function -{S}(x::Poly{ZZ, S}, y::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_sub_si, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

function -{S}(x::Poly{ZZ, S}, y::ZZ)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_sub_fmpz, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{ZZ}), 
               &(z.data), &(x.data), &y)
   return z
end

function -{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Int)
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_sub_si, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

function -{S, M}(x::Poly{Residue{ZZ, M}, S}, y::ZZ)
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_sub_fmpz, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{ZZ}), 
               &(z.data), &(x.data), &y)
   return z
end

-{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Residue{ZZ, M}) = x - y.data

function -{S}(x::Int, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_si_sub, :libflint), Void, 
                (Ptr{fmpz_poly}, Int, Ptr{fmpz_poly}), 
               &(z.data), x, &(y.data))
   return z
end

function -{S, M}(x::Int, y::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_si_sub, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Int, Ptr{fmpz_mod_poly}), 
               &(z.data), x, &(y.data))
   return z
end

function -{S}(x::ZZ, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_fmpz_sub, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{ZZ}, Ptr{fmpz_poly}), 
               &(z.data), &x, &(y.data))
   return z
end

function -{S, M}(x::ZZ, y::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_fmpz_sub, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{ZZ}, Ptr{fmpz_mod_poly}), 
               &(z.data), &x, &(y.data))
   return z
end

-{S, M}(x::Residue{ZZ, M}, y::Poly{Residue{ZZ, M}, S}) = x.data - y

+{S}(x::Int, y::Poly{ZZ, S}) = y + x

+{S}(x::ZZ, y::Poly{ZZ, S}) = y + x

*{S}(x::Poly{ZZ, S}, y::Int) = y*x

*{S}(x::Poly{ZZ, S}, y::ZZ) = y*x

*{T <: Ring, S}(a::Poly{T, S}, b::Int) = b*a

*{T <: Ring, S}(a::Poly{T, S}, b::ZZ) = b*a

###########################################################################################
#
#   Truncation
#
###########################################################################################

function truncate{S}(a::Poly{ZZ, S}, n::Int)
   n < 0 && throw(DomainError())
   
   if a.data.length <= n
      return a
   end

   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_set_trunc, :libflint), Void,
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int),
               &(z.data), &(a.data), n)

   return z
end

function truncate{S, M}(a::Poly{Residue{ZZ, M}, S}, n::Int)
   n < 0 && throw(DomainError())
   
   if a.data.length <= n
      return a
   end

   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_set_trunc, :libflint), Void,
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Int),
               &(z.data), &(a.data), n)
   return z
end

function truncate{T <: Ring, S}(a::Poly{T, S}, n::Int)
   n < 0 && throw(DomainError())
   
   lena = a.data.length

   if lena <= n
      return a
   end

   lenz = min(lena, n)
   d = Array(T, lenz)

   for i = 1:lenz
      d[i] = a.data.coeffs[i]
   end

   z = Poly(Poly{T, S}, d)

   z.data.length = normalise(z, lenz)

   return z
end

function mullow{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S}, n::Int)
   n < 0 && throw(DomainError())
   
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_mullow, :libflint), Void,
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int),
               &(z.data), &(x.data), &(y.data), n)
   return z
end

function mullow{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S}, n::Int)
   n < 0 && throw(DomainError())
   
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_mullow, :libflint), Void,
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Int),
               &(z.data), &(x.data), &(y.data), n)
   return z
end

function mullow{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S}, n::Int)
   lena = a.data.length
   lenb = b.data.length

   if lena == 0 || lenb == 0
      return zero(Poly{T, S})
   end

   if n < 0
      n = 0
   end

   t = T()

   lenz = min(lena + lenb - 1, n)

   d = Array(T, lenz)

   for i = 1:min(lena, lenz)
      d[i] = a.data.coeffs[i]*b.data.coeffs[1]
   end

   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = a.data.coeffs[lena]*b.data.coeffs[j]
      end
   end

   z = Poly(Poly{T, S}, d)

   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            mul!(t, a.data.coeffs[i], b.data.coeffs[j])
            addeq!(z.data.coeffs[i + j - 1], t)
         end
      end
   end
        
   z.data.length = normalise(z, lenz)

   return z
end

###########################################################################################
#
#   Reversal
#
###########################################################################################

function reverse{T <: Ring, S}(x::Poly{T, S}, len::Int)
   len < 0 && throw(DomainError())
   v = Array(T, len)
   for i = 1:len
      v[i] = coeff(x, len - i)
   end
   r = Poly(Poly{T, S}, v)
   r.data.length = normalise(r, len)
   return r
end

function reverse{S}(x::Poly{ZZ, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_reverse, :libflint), Void,
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int),
               &(z.data), &(x.data), len)
   return z
end

function reverse{S, M}(x::Poly{Residue{ZZ, M}, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_reverse, :libflint), Void,
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Int),
               &(z.data), &(x.data), len)
   return z
end

function reverse{T <: Ring, S}(x::Poly{T, S})
   reverse(x, x.data.length)
end

###########################################################################################
#
#   Shifting
#
###########################################################################################

function shift_left{T <: Ring, S}(x::Poly{T, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = x.data.length
   v = Array(T, xlen + len)
   for i = 1:len
      v[i] = zero(T)
   end
   for i = 1:xlen
      v[i + len] = coeff(x, i - 1)
   end
   return Poly(Poly{T, S}, v)
end

function shift_right{T <: Ring, S}(x::Poly{T, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = x.data.length
   if len >= xlen
      return zero(Poly{T, S})
   end
   v = Array(T, xlen - len)
   for i = 1:xlen - len
      v[i] = coeff(x, i + len - 1)
   end
   return Poly(Poly{T, S}, v)
end

function shift_left{S, M}(x::Poly{Residue{ZZ, M}, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_shift_left, :libflint), Void,
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Int),
               &(z.data), &(x.data), len)
   return z
end

function shift_right{S, M}(x::Poly{Residue{ZZ, M}, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_shift_right, :libflint), Void,
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Int),
               &(z.data), &(x.data), len)
   return z
end

function shift_left{S}(x::Poly{ZZ, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_shift_left, :libflint), Void,
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int),
               &(z.data), &(x.data), len)
   return z
end

function shift_right{S}(x::Poly{ZZ, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_shift_right, :libflint), Void,
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int),
               &(z.data), &(x.data), len)
   return z
end


###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{S}(x::Poly{ZZ, S}, y::Int)
   y < 0 && throw(DomainError())
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_pow, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

function ^{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Int)
   y < 0 && throw(DomainError())
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_pow, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

function ^{T <: Ring, S}(a::Poly{T, S}, b::Int)
   b < 0 && throw(DomainError())
   # special case powers of x for constructing polynomials efficiently
   if a.data.length == 2 && a.data.coeffs[1] == 0 && a.data.coeffs[2] == 1
      d = Array(T, b + 1)
      d[b + 1] = a.data.coeffs[2]
      for i = 1:b
         d[i] = a.data.coeffs[1]
      end
      z = Poly(Poly{T, S}, d)
      z.data.length = b + 1
      return z
   elseif a.data.length == 0
      return zero(Poly{T, S})
   elseif a.data.length == 1
      return Poly(Poly{T, S}, [a.data.coeffs[1]^b])
   elseif b == 0
      return one(Poly{T, S})
   else
      bit = ~((~uint(0)) >> 1)
      while (int(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = z*z
         if (int(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end
   
###########################################################################################
#
#   Modular arithmetic
#
###########################################################################################

function mulmod{T <: Union(Residue, Field), S}(a::Poly{T, S}, b::Poly{T, S}, d::Poly{T, S})
   return mod(a*b, d)
end

function mulmod{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S}, f::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_mulmod, :libflint), Void,
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}),
               &(z.data), &(x.data), &(y.data), &(f.data))
   return z
end

function powmod{T <: Union(Residue, Field), S}(a::Poly{T, S}, b::Int, d::Poly{T, S})
   if a.data.length == 0
      return zero(Poly{T, S})
   elseif a.data.length == 1
      return Poly(Poly{T, S}, [a.data.coeffs[1]^b])
   elseif b == 0
      return one(Poly{T, S})
   else
      if b < 0
         a = invmod(a, d)
         b = -b
      end
      bit = ~((~uint(0)) >> 1)
      while (int(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = mulmod(z, z, d)
         if (int(bit) & b) != 0
            z = mulmod(z, a, d)
         end
         bit >>= 1
      end
      return z
   end
end

function powmod{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Int, f::Poly{Residue{ZZ, M}, S})
   if y < 0
      x = invmod(x, f)
      y = -y
   end
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_powmod_ui_binexp, :libflint), Void,
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Int, Ptr{fmpz_mod_poly}),
               &(z.data), &(x.data), y, &(f.data))
   return z
end

function invmod{T <: Union(Residue, Field), S}(a::Poly{T, S}, b::Poly{T, S})
   g, z = gcdinv(a, b)
   if g != 1
      error("Impossible inverse in invmod")
   end
   return z
end

function invmod{S, M}(x::Poly{Residue{ZZ, M}, S}, f::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_invmod, :libflint), Void,
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}),
               &(z.data), &(x.data), &(f.data))
   return z
end


###########################################################################################
#
#   Comparisons
#
###########################################################################################

function =={S}(x::Poly{ZZ, S}, y::ZZ) 
   if x.data.length > 1
      return false
   elseif x.data.length == 1 
      z = ZZ();
      ccall((:fmpz_poly_get_coeff_fmpz, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_poly}, Int), 
               &z, &(x.data), 0)
      return z == y
   else
      return y == 0
   end 
end

=={S}(x::Poly{ZZ, S}, y::Int) = x == ZZ(y)

function =={S, M}(x::Poly{Residue{ZZ, M}, S}, y::ZZ) 
   if x.data.length > 1
      return false
   elseif x.data.length == 1 
      z = ZZ();
      ccall((:fmpz_poly_get_coeff_fmpz, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_mod_poly}, Int), 
               &z, &(x.data), 0)
      return z == y
   else
      return y == 0
   end 
end

=={S, M}(x::Poly{Residue{ZZ, M}, S}, y::Int) = x == ZZ(y)

=={S, M}(x::Poly{Residue{ZZ, M}, S}, y::Residue{ZZ, M}) = x == y.data

=={T <: Ring, S}(x::Poly{T, S}, y::Int) = ((x.data.length == 0 && y == 0)
                                        || (x.data.length == 1 && x.data.coeffs[1] == y))

=={T <: Ring, S}(x::Poly{T, S}, y::ZZ) = ((x.data.length == 0 && y == 0)
                                        || (x.data.length == 1 && x.data.coeffs[1] == y))

=={S}(x::Poly{ZZ, S}, y::Poly{ZZ, S}) = ccall((:fmpz_poly_equal, :libflint), Bool, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}), &(x.data), &(y.data))

=={S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S}) = ccall((:fmpz_mod_poly_equal, :libflint), Bool, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), &(x.data), &(y.data))

function =={T<: Ring, S}(x::Poly{T, S}, y::Poly{T, S})
   if x.data.length != y.data.length
      return false
   else
      for i = 1:x.data.length
         if x.data.coeffs[i] != y.data.coeffs[i]
            return false
         end
      end
   end
   return true
end

=={T<: Ring, S}(x::Int, y::Poly{T, S}) = y == x

=={T<: Ring, S}(x::ZZ, y::Poly{T, S}) = y == x

=={S, M}(x::Residue{ZZ, M}, y::Poly{Residue{ZZ, M}, S}) = y == x

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{T <: Ring, S}(a::Poly{T, S}, b::T)
   b == 0 && throw(DivideError())
   d = Array(T, a.data.length)
   for i = 1:a.data.length
      d[i] = divexact(a.data.coeffs[i], b)
   end
   z = Poly(Poly{T, S}, d)
   z.data.length = a.data.length
   return z
end

function divexact{T <: Ring, S}(a::Poly{T, S}, b::ZZ)
   b == 0 && throw(DivideError())
   d = Array(T, a.data.length)
   for i = 1:a.data.length
      d[i] = divexact(a.data.coeffs[i], b)
   end
   z = Poly(Poly{T, S}, d)
   z.data.length = a.data.length
   return z
end

function divexact{T <: Ring, S}(a::Poly{T, S}, b::Int)
   b == 0 && throw(DivideError())
   d = Array(T, a.data.length)
   for i = 1:a.data.length
      d[i] = divexact(a.data.coeffs[i], b)
   end
   z = Poly(Poly{T, S}, d)
   z.data.length = a.data.length
   return z
end

function divexact{S}(x::Poly{ZZ, S}, y::ZZ)
   y == 0 && throw(DivideError())
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_scalar_divexact_fmpz, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{ZZ}), 
               &(z.data),  &(x.data), &y)
   return z
end

function divexact{S}(x::Poly{ZZ, S}, y::Int)
   y == 0 && throw(DivideError())
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_scalar_divexact_si, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), 
               &(z.data),  &(x.data), y)
   return z
end

function divexact{S, M}(x::Poly{Residue{ZZ, M}, S}, y::ZZ)
   y == 0 && throw(DivideError())
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_scalar_div_fmpz, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{ZZ}), 
               &(z.data),  &(x.data), &y)
   return z
end

divexact{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Int) = divexact(x, ZZ(y))

divexact{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Residue{ZZ, M}) = divexact(x, y.data)

function divexact{T <: Ring, S}(f::Poly{T, S}, g::Poly{T, S})
   g == 0 && throw(DivideError())
   if f == 0
      return zero(Poly{T, S})
   end
   lenq = f.data.length - g.data.length + 1
   d = Array(T, lenq)
   for i = 1:lenq
      d[i] = zero(T)
   end
   q = Poly(Poly{T, S}, d)
   x = gen(Poly{T, S})
   leng = g.data.length
   while f.data.length >= leng
      lenf = f.data.length
      q1 = q.data.coeffs[lenf - leng + 1] = divexact(f.data.coeffs[lenf], g.data.coeffs[leng])
      f = f - q1*g*x^(lenf - leng)
   end
   q.data.length = lenq
   return q
end

function divexact{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   y == 0 && throw(DivideError())
   if x == 0
      return zero(Poly{ZZ, S})
   end
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_div, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data),  &(x.data), &(y.data))
   return z
end

function divexact{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   y == 0 && throw(DivideError())
   if x == 0
      return zero(Poly{Residue{ZZ, M}, S})
   end
   q = Poly{Residue{ZZ, M}, S}()
   r = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_divrem, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(q.data), &(r.data), &(x.data), &(y.data))
   return q
end

###########################################################################################
#
#   Euclidean division
#
###########################################################################################

function mod{T <: Union(Field, Residue), S}(f::Poly{T, S}, g::Poly{T, S})
   if g.data.length == 0
      raise(DivideError())
   end
   if f.data.length >= g.data.length
      b = g.data.coeffs[g.data.length]
      g = inv(b)*g
      x = gen(Poly{T, S})
      while f.data.length >= g.data.length
         f -= coeff(f, f.data.length - 1)*g*x^(f.data.length - g.data.length)
      end
   end
   return f
end

function mod{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   y == 0 && throw(DivideError())
   r = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_rem, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(r.data), &(x.data), &(y.data))
   return r
end

function divrem{T <: Union(Field, Residue), S}(f::Poly{T, S}, g::Poly{T, S})
   if g.data.length == 0
      raise(DivideError())
   end
   if f.data.length < g.data.length
      return zero(Poly{T, S}), f
   end
   binv = inv(lead(g))
   g = binv*g
   x = gen(Poly{T, S})
   qlen = f.data.length - g.data.length + 1
   d = Array(T, qlen)
   for i = 1:qlen
      d[i] = zero(T)
   end
   q = Poly(Poly{T, S}, d)
   while f.data.length >= g.data.length
      q1 = coeff(f, f.data.length - 1)
      setcoeff!(q, f.data.length - g.data.length, q1*binv)
      f -= q1*g*x^(f.data.length - g.data.length)
   end
   return q, f
end

function divrem{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   y == 0 && throw(DivideError())
   if x == 0
      return zero(Poly{Residue{ZZ, M}, S})
   end
   q = Poly{Residue{ZZ, M}, S}()
   r = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_divrem, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(q.data), &(r.data), &(x.data), &(y.data))
   return q, r
end

###########################################################################################
#
#   Pseudodivision
#
###########################################################################################

function pseudorem{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   y == 0 && throw(DivideError())
   z = Poly{ZZ, S}()
   d = 0
   ccall((:fmpz_poly_pseudo_rem, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{Int}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &d, &(x.data), &(y.data))
   return z
end

function pseudodivrem{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   y == 0 && throw(DivideError())
   q = Poly{ZZ, S}()
   r = Poly{ZZ, S}()
   d = 0
   ccall((:fmpz_poly_pseudo_divrem_divconquer, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{Int}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(q.data), &(r.data), &d, &(x.data), &(y.data))
   return (q, r)
end

function pseudorem{T <: Ring, S}(f::Poly{T, S}, g::Poly{T, S})
   g == 0 && throw(DivideError())
   b = coeff(g, g.data.length - 1)
   x = gen(Poly{T, S})
   while f.data.length >= g.data.length
      f = f*b - coeff(f, f.data.length - 1)*g*x^(f.data.length - g.data.length)
   end
   return f
end

function pseudodivrem{T <: Ring, S}(f::Poly{T, S}, g::Poly{T, S})
   g == 0 && throw(DivideError())
   if f.data.length < g.data.length
      return zero(Poly{T, S}), f
   end
   lenq = f.data.length - g.data.length + 1
   v = Array(T, lenq)
   for i = 1:lenq
      v[i] = zero(T)
   end
   q = Poly(Poly{T, S}, v)
   b = coeff(g, g.data.length - 1)
   x = gen(Poly{T, S})
   while f.data.length >= g.data.length
      for i = f.data.length - g.data.length + 2:lenq
         setcoeff!(q, i - 1, coeff(q, i - 1) * b)
      end
      setcoeff!(q, f.data.length - g.data.length, coeff(f, f.data.length - 1))
      f = f*b - coeff(f, f.data.length - 1)*g*x^(f.data.length - g.data.length)
   end
   q.data.length = normalise(q, lenq)
   return q, f
end

###########################################################################################
#
#   Content, primitive part, GCD and LCM
#
###########################################################################################

function gcd{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   if b.data.length > a.data.length
      (a, b) = (b, a)
   end
   if b == 0
      return a
   end
   c = gcd(content(a), content(b))
   a = divexact(a, c)
   b = divexact(b, c)
   g = one(Poly{T, S})
   h = one(Poly{T, S})
   while true
      d = a.data.length - b.data.length
      r = pseudorem(a, b)
      if r == 0
         break
      end
      if r.data.length == 1
         b = one(Poly{T, S})
         break
      end
      (a, b) = (b, divexact(r, g*h^d))
      g = lead(a)
      if d > 1
         h = divexact(g^d, h^(d - 1))
      else
         h = h^(1 - d)*g^d
      end
   end
   return c*primpart(b)
end

function gcd{T <: Union(Field, Residue), S}(a::Poly{T, S}, b::Poly{T, S})
   if a.data.length > b.data.length
      (a, b) = (b, a)
   end
   if b == 0
      return a
   end
   g = gcd(content(a), content(b))
   a = divexact(a, g)
   b = divexact(b, g)
   while a != 0
      (a, b) = (mod(b, a), a)
   end
   b = g*b
   return inv(lead(b))*b
end

function gcd{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_gcd, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function gcd{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_gcd, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function lcm{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   return a*divexact(b, gcd(a, b))
end

function content{T <: Ring, S}(a::Poly{T, S})
   z = coeff(a, 0)
   for i = 2:a.data.length
      z = gcd(z, coeff(a, i - 1))
   end
   return z
end

function content{S}(x::Poly{ZZ, S})
   z = ZZ()
   ccall((:fmpz_poly_content, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_poly}), 
               &z, &(x.data))
   return z
end

function primpart{T <: Ring, S}(a::Poly{T, S})
   d = content(a)
   return divexact(a, d)
end

function primpart{S}(x::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_primitive_part, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data))
   return z
end

###########################################################################################
#
#   Evaluation/composition
#
###########################################################################################

function evaluate{S}(x::Poly{ZZ, S}, y::ZZ)
   z = ZZ()
   ccall((:fmpz_poly_evaluate_fmpz, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_poly}, Ptr{ZZ}), 
               &z, &(x.data), &y)
   return z
end

function evaluate{S}(x::Poly{ZZ, S}, y::Int)
   z = ZZ()
   ccall((:fmpz_poly_evaluate_fmpz, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_poly}, Ptr{ZZ}), 
               &z, &(x.data), &ZZ(y))
   return z
end

function evaluate{S, M}(x::Poly{Residue{ZZ, M}, S}, y::ZZ)
   z = ZZ()
   ccall((:fmpz_mod_poly_evaluate_fmpz, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_mod_poly}, Ptr{ZZ}), 
               &z, &(x.data), &y)
   return z
end

evaluate{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Residue{ZZ, M}) = evaluate(x, y.data)

function evaluate{T <: Ring, S}(a::Poly{T, S}, b::T)
   i = a.data.length
   if i == 0
       return zero(T)
   end
   z = a.data.coeffs[i]
   while i > 1
      i -= 1
      z = z*b + a.data.coeffs[i]
   end
   return z
end

function evaluate{T <: Ring, S}(a::Poly{T, S}, b::Int)
   return evaluate(a, convert(T, b))
end

function evaluate{T <: Ring, S}(a::Poly{T, S}, b::ZZ)
   return evaluate(a, convert(T, b))
end

function compose{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_compose, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function compose{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_compose, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function compose{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   i = a.data.length
   if i == 0
       return zero(T)
   end
   z = a.data.coeffs[i]
   while i > 1
      i -= 1
      z = z*b + a.data.coeffs[i]
   end
   return z
end

###########################################################################################
#
#   Derivative
#
###########################################################################################

function deriv{S}(x::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_derivative, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data))
   return z
end

function deriv{S, M}(x::Poly{Residue{ZZ, M}, S})
   z = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_derivative, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(z.data), &(x.data))
   return z
end

function deriv{T <: Ring, S}(a::Poly{T, S})
   if a == 0
      return zero(Poly{T, S})
   end
   len = a.data.length
   d = Array(T, len - 1)
   for i = 1:len - 1
      d[i] = i*a.data.coeffs[i + 1]
   end
   z = Poly(Poly{T, S}, d)
   z.data.length = normalise(z, len - 1)
   return z
end

###########################################################################################
#
#   Integral
#
###########################################################################################

function integral{T <: Union(Field, Residue), S}(x::Poly{T, S})
   len = x.data.length
   v = Array(T, len + 1)
   v[1] = zero(T)
   for i = 1:len
      v[i + 1] = divexact(coeff(x, i - 1), T(i))
   end
   p = Poly(Poly{T, S}, v)
   len = len + 1
   while len > 0 && coeff(p, len - 1) == 0 # cannot use normalise here
      len -= 1
   end
   return p
end

###########################################################################################
#
#   Resultant
#
###########################################################################################

function resultant{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   z = ZZ()
   ccall((:fmpz_poly_resultant, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &z, &(x.data), &(y.data))
   return z
end

function resultant{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   if a.data.length == 0 || b.data.length == 0
      return zero(T)
   end
   sgn = 1
   if a.data.length < b.data.length
      a, b = b, a
      if iseven(a.data.length) && iseven(b.data.length)
         sgn = -sgn
      end
   end
   lena = a.data.length
   lenb = b.data.length
   if lenb == 1
      return b.data.coeffs[1]^(lena - 1)
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   g = one(T)
   h = one(T)
   while lenb > 1
      d = lena - lenb
      if iseven(lena) && iseven(lenb)
         sgn = -sgn
      end
      B, A = pseudorem(A, B), B
      lena = lenb
      lenb = B.data.length
      if lenb == 0
         return zero(T) 
      end
      s = h^d
      B = divexact(B, g*s)
      g = lead(A)
      h = divexact(h*g^d, s)
   end
   s = divexact(h*lead(B)^(lena - 1), h^(lena - 1))
   res = c1^(lenb - 1)*c2^(lena - 1)*s*sgn
end

function resultant{T <: Residue, S}(a::Poly{T, S}, b::Poly{T, S})
   if a.data.length == 0 || b.data.length == 0
      return zero(T)
   end
   sgn = 1
   if a.data.length < b.data.length
      a, b = b, a
      if iseven(a.data.length) && iseven(b.data.length)
         sgn = -sgn
      end
   end
   lena = a.data.length
   lenb = b.data.length
   if lenb == 1
      return b.data.coeffs[1]^(lena - 1)
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   s = 1
   while lenb > 1
      if iseven(lena) && iseven(lenb)
         sgn = -sgn
      end
      B, A = mod(A, B), B
      s *= lead(A)^(lena - B.data.length)
      lena = lenb
      lenb = B.data.length
      if lenb == 0
         return zero(T) 
      end
   end
   s *= lead(B)^(lena - 1)
   res = c1^(lenb - 1)*c2^(lena - 1)*s*sgn
end

function resultant{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   z = ZZ()
   ccall((:fmpz_mod_poly_resultant, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &z, &(x.data), &(y.data))
   return Residue{ZZ, M}(z)
end

###########################################################################################
#
#   Discriminant
#
###########################################################################################

function discriminant{S}(x::Poly{ZZ, S})
   z = ZZ()
   ccall((:fmpz_poly_discriminant, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_poly}), 
               &z, &(x.data))
   return z
end

function discriminant{S, M}(x::Poly{Residue{ZZ, M}, S})
   z = ZZ()
   ccall((:fmpz_mod_poly_discriminant, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_mod_poly}), 
               &z, &(x.data))
   return Residue{ZZ, M}(z)
end

function discriminant{T <: Ring, S}(a::Poly{T, S})
   d = deriv(a)
   z = resultant(a, d)
   if a.data.length - d.data.length == 1
      z = divexact(z, lead(a))
   else
      z = z*lead(a)^(a.data.length - d.data.length - 2)
   end
   mod4 = (a.data.length + 3)%4 # degree mod 4
   return mod4 == 2 || mod4 == 3 ? -z : z
end

###########################################################################################
#
#   Bezout
#
###########################################################################################

function bezout{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   z = ZZ()
   u = Poly{ZZ, S}()
   v = Poly{ZZ, S}()
   ccall((:fmpz_poly_xgcd_modular, :libflint), Void, 
                (Ptr{ZZ}, Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &z, &(u.data), &(v.data), &(x.data), &(y.data))
   return (z, u, v)
end

function bezout{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   g = Poly{Residue{ZZ, M}, S}()
   s = Poly{Residue{ZZ, M}, S}()
   t = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_xgcd, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(g.data), &(s.data), &(t.data), &(x.data), &(y.data))
   return (g, s, t)
end

function bezout{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   if a.data.length == 0 || b.data.length == 0
      return zero(T), zero(Poly{T, S}), zero(Poly{T, S})
   end
   sgn = 1
   swap = false
   if a.data.length < b.data.length
      a, b = b, a
      swap = true
      if iseven(a.data.length) && iseven(b.data.length)
         sgn = -sgn
      end
   end
   lena = a.data.length
   lenb = b.data.length
   if lenb == 1
      s1 = zero(T)
      t1 = one(T)
      r1 = b.data.coeffs[1]^(lena - 1)
      if swap
         s1, t1 = t1, s1
      end
      if sgn
         s1, t1, r1 = -s1, -t1, -r1
      end
      return r1, convert(Poly{T, S}, s1), convert(Poly{T, S}, t1)
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   g = one(T)
   h = one(T)
   u1, u2 = one(Poly{T, S}), zero(Poly{T, S})
   v1, v2 = zero(Poly{T, S}), one(Poly{T, S})
   while lenb > 1
      d = lena - lenb
      if iseven(lena) && iseven(lenb)
         sgn = -sgn
      end
      (Q, B), A = pseudodivrem(A, B), B
      lena = lenb
      lenb = B.data.length
      if lenb == 0
         return zero(T), zero(Poly{T, S}), zero(Poly{T, S})
      end
      s = h^d
      B = divexact(B, g*s)
      t = lead(A)^(d + 1)
      u2, u1 = divexact(u1*t - Q*u2, g*s), u2
      v2, v1 = divexact(v1*t - Q*v2, g*s), v2 
      g = lead(A)
      h = divexact(h*g^d, s)
   end
   s = divexact(h*lead(B)^(lena - 1), h^(lena - 1))
   res = c1^(lenb - 1)*c2^(lena - 1)*s*sgn
   if swap
      u2, v2 = v2, u2
   end
   return res, u2, v2
end

function bezout{T <: Union(Field, Residue), S}(a::Poly{T, S}, b::Poly{T, S})
   if a.data.length == 0
      return b, zero(Poly{T, S}), one(Poly{T, S})
   end
   if b.data.length == 0
      return a, one(Poly{T, S}), zero(Poly{T, S})
   end
   swap = false
   if a.data.length < b.data.length
      a, b = b, a
      swap = true
   end
   lena = a.data.length
   lenb = b.data.length
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   u1, u2 = inv(c1), zero(Poly{T, S})
   v1, v2 = zero(Poly{T, S}), inv(c2)
   while lenb > 0
      d = lena - lenb
      (Q, B), A = divrem(A, B), B
      lena = lenb
      lenb = B.data.length
      u2, u1 = u1 - Q*u2, u2
      v2, v1 = v1 - Q*v2, v2 
   end
   if swap
      u1, v1 = v1, u1
   end
   d = gcd(c1, c2)
   A, u1, v1 = d*A, d*u1, d*v1
   d = inv(lead(A))
   return d*A, d*u1, d*v1
end

function gcdinv{T <: Union(Field, Residue), S}(a::Poly{T, S}, b::Poly{T, S})
   if a.data.length == 0
      if b.data.length == 0
         return 0, zero(Poly{T, S})
      else
         d = inv(lead(b))
         return b*d, zero(Poly{T, S})
      end
   end
   if b.data.length == 0
      d = inv(lead(b))
      return a*d, d
   end
   if a.data.length < b.data.length
      a, b = b, a
      u1, u2 = zero(Poly{T, S}), one(Poly{T, S})
   else
      u1, u2 = one(Poly{T, S}), zero(Poly{T, S})
   end
   lena = a.data.length
   lenb = b.data.length
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   u1 *= inv(c1)
   u2 *= inv(c2)
   while lenb > 0
      d = lena - lenb
      (Q, B), A = divrem(A, B), B
      lena = lenb
      lenb = B.data.length
      u2, u1 = u1 - Q*u2, u2
   end
   d = gcd(c1, c2)
   A, u1 = d*A, d*u1
   d = inv(lead(A))
   return d*A, d*u1
end

function gcdinv{S, M}(x::Poly{Residue{ZZ, M}, S}, y::Poly{Residue{ZZ, M}, S})
   g = Poly{Residue{ZZ, M}, S}()
   s = Poly{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_gcdinv, :libflint), Void, 
                (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), 
               &(g.data), &(s.data), &(x.data), &(y.data))
   return (g, s)
end

###########################################################################################
#
#   Special polynomials
#
###########################################################################################

function chebyshev_t{S}(::Type{Poly{ZZ, S}}, n::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_chebyshev_t, :libflint), Void, (Ptr{fmpz_poly}, Int), &(z.data), n)
   return z
end
   
function chebyshev_u{S}(::Type{Poly{ZZ, S}}, n::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_chebyshev_u, :libflint), Void, (Ptr{fmpz_poly}, Int), &(z.data), n)
   return z
end

function chebyshev_t_pair{S <: Ring}(n::Int, x::S)
   if n == 0
      return one(S), x
   elseif n == 1
      return x, one(S)
   elseif n < 0
      a, b = chebyshev_t_pair(1-n, x)
      return b, a
   elseif iseven(n)
      a, b = chebyshev_t_pair(n>>1, x)
      return 2*(a*a) - 1, 2*(a*b) - x
   else
      a, b = chebyshev_t_pair((n>>1)+1, x)
      return 2*(a*b) - x, 2*(b*b) - 1
   end
end

function chebyshev_t{S <: Ring}(n::Int, x::S)
   if n == 0
      return one(S)
   elseif n == 1
      return x
   elseif n < 0
      return chebyshev_t(-n, x)
   elseif iseven(n)
      a = chebyshev_t(n>>1, x)
      return 2*(a*a) - 1
   else
      a, b = chebyshev_t_pair((n>>1)+1, x)
      return 2*(a*b) - x
   end
end

function chebyshev_u_pair{S <: Ring}(n::Int, x::S)
   if n == 0
      return one(S), zero(S)
   elseif n == 1
      return 2*x, one(S)
   elseif n == -1
      return zero(S), -one(S)
   elseif n < -1
      a, b = chebyshev_u_pair(-1-n, x)
      return -b, -a
   elseif iseven(n)
      a, b = chebyshev_u_pair(n>>1, x)
      return (a+b)*(a-b), 2*b*(a-x*b)
   else
      a, b = chebyshev_u_pair(n>>1, x)
      return 2*a*(x*a-b), (a+b)*(a-b)
   end
end

function chebyshev_u{S <: Ring}(n::Int, x::S)
   if n == 0
      return one(S)
   elseif n == 1
      return 2*x
   elseif n == -1
      return zero(S)
   elseif n < -1
      return -chebyshev_u(-2-n, x)
   elseif iseven(n)
      a, b = chebyshev_u_pair(n>>1, x)
      return (a+b)*(a-b)
   else
      a, b = chebyshev_u_pair(n>>1, x)
      return 2*a*(x*a-b)
   end
end

function cyclotomic{S}(::Type{Poly{ZZ, S}}, n::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_cyclotomic, :libflint), Void, (Ptr{fmpz_poly}, Int), &(z.data), n)
   return z
end
   
function swinnerton_dyer{S}(::Type{Poly{ZZ, S}}, n::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_swinnerton_dyer, :libflint), Void, 
                                                    (Ptr{fmpz_poly}, Int), &(z.data), n)
   return z
end
   
function cos_minpoly{S}(::Type{Poly{ZZ, S}}, n::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_cos_minpoly, :libflint), Void, (Ptr{fmpz_poly}, Int), &(z.data), n)
   return z
end
   
function theta_qexp{S}(::Type{Poly{ZZ, S}}, e::Int, n::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_theta_qexp, :libflint), Void, 
                                            (Ptr{fmpz_poly}, Int, Int), &(z.data), e, n)
   return z
end

function eta_qexp{S}(::Type{Poly{ZZ, S}}, e::Int, n::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_eta_qexp, :libflint), Void, 
                                            (Ptr{fmpz_poly}, Int, Int), &(z.data), e, n)
   return z
end
   
###########################################################################################
#
#   PolynomialRing constructor
#
###########################################################################################

function PolynomialRing{T <: Ring}(::Type{T}, s::String)
   S = symbol(s)
   T1 = Poly{T, S}
   T2 = T
   
   # Conversions and promotions

   Base.convert(::Type{T1}, x::T) = Poly(T1, [x])
   Base.promote_rule(::Type{T1}, ::Type{T}) = T1

   P = T2.parameters
   while length(P) > 0
      T2 = P[1]
      if isa(T2, DataType) && T2 <: Ring
         Base.convert(::Type{T1}, x::T2) = Poly(T1, [convert(T, x)])
         Base.promote_rule(::Type{T1}, ::Type{T2}) = T1
         P = T2.parameters
      else
         break
      end
   end

   Base.convert(::Type{T1}, x::Integer) = Poly(T1, [convert(T, x)])
   Base.promote_rule{R <: Integer}(::Type{T1}, ::Type{R}) = T1

   # (Type, gen) 

   return (Poly{T, S}, Poly(Poly{T, S}, [T(0), T(1)]))
end
