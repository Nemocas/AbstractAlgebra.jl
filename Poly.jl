export Poly, PolynomialRing, coeff, zero, one, gen, isgen, normalise, chebyshev_t,
       chebyshev_u, theta_qexp, eta_qexp, swinnerton_dyer, cos_minpoly, cyclotomic,
       pseudorem, pseudodivrem, primpart, content, divexact, evaluate, compose, deriv,
       resultant, lead, discriminant, bezout, truncate, mullow

import Base: convert, zero

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

type fmpz_poly <: Ring
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
end

type PolyStruct{T <: Ring}
   coeffs :: Array{T, 1}
   length :: Int
end

type Poly{T <: Ring, S} <: Ring
   data :: Union(fmpz_poly, PolyStruct{T})
   
   function Poly(a :: Array{T, 1})
      if T == ZZ
         d = new(fmpz_poly(C_NULL, 0, 0))
         ccall((:fmpz_poly_init2, :libflint), Void, (Ptr{fmpz_poly}, Int), &(d.data), length(a))
         finalizer(d, _fmpz_poly_clear_fn)
         for i = 1:length(a)
            ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Void, (Ptr{fmpz_poly}, Int, Ptr{ZZ}),
                 &(d.data), i - 1, &a[i])
         end
      else
         len = length(a)
         while len > 0 && a[len] == 0
            len -= 1
         end
         d = new(PolyStruct(a, len))
      end
      return d
   end   

   Poly() = Poly{T, S}(Array(T, 0))
   Poly(a::Integer) = a == 0 ? Poly{T, S}(Array(T, 0)) : Poly{T, S}([T(a)])
   Poly(a::T) = Poly{T, S}([a])
   Poly(a::Poly{T, S}) = a
   Poly{R <: Ring}(a::R) = convert(Poly{T, S}, a)
end

function _fmpz_poly_clear_fn(a :: Poly{ZZ})
   ccall((:fmpz_poly_clear, :libflint), Void, (Ptr{fmpz_poly},), &(a.data))
end
   
###########################################################################################
#
#   Basic manipulation
#
###########################################################################################    
   
function normalise{T <: Ring, S}(a::Poly{T, S}, len::Int)
   while len > 0 && a.data.coeffs[len] == 0
      len -= 1
   end

   return len
end

# Julia wants "objects" to be immutable, so set_coeff is not provided

function coeff{S}(x::Poly{ZZ, S}, n::Int)
   z = ZZ()
   ccall((:fmpz_poly_get_coeff_fmpz, :libflint), Void, (Ptr{ZZ}, Ptr{fmpz_poly}, Int), &z, &(x.data), n)
   return z
end

coeff{T <: Ring, S}(a::Poly{T, S}, n::Int) = n >= a.data.length ? 0 : a.data.coeffs[n + 1]

lead{S}(x::Poly{ZZ, S}) = x.data.length == 0 ? zero(ZZ) : coeff(x, x.data.length - 1)

lead{T <: Ring, S}(a::Poly{T, S}) = a.data.length == 0 ? zero(T) : a.data.coeffs[a.data.length]

isgen{S}(x::Poly{ZZ, S}) = bool(ccall((:fmpz_poly_is_x, :libflint), Int, (Ptr{fmpz_poly},), &(x.data)))

isgen{T <: Ring, S}(a::Poly{T, S}) = a.data.length == 2 && a.data.coeffs[1] == 0 && a.data.coeffs[2] == 1

zero{S}(::Type{Poly{ZZ, S}}) = Poly{ZZ, S}(0)

zero{T <: Ring, S}(::Type{Poly{T, S}}) = Poly{T, S}(0)

one{S}(::Type{Poly{ZZ, S}}) = Poly{ZZ, S}(1)

one{T <: Ring, S}(::Type{Poly{T, S}}) = Poly{T, S}(1)

gen{S}(::Type{Poly{ZZ, S}}) = Poly{ZZ, S}([ZZ(0), 1])

gen{T <: Ring, S}(::Type{Poly{T, S}}) = Poly{T, S}([T(0), T(1)])

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

function show{T <: Ring, S}(io::IO, x::Poly{T, S})
   len = x.data.length

   if len == 0
      print(io, zero(T))
   else
      for i = 1:len - 1
         c = x.data.coeffs[len - i + 1]
         bracket = isa(c, Poly) && c.data.length != 1
         if c != 0
            if i != 1
               print(io, "+")
            end
            if c != 1
               if bracket
                  print(io, "(")
               end
               show(io, c)
               if bracket
                  print(io, ")")
               end
               print(io, "*")
            end
            print(io, string(S))
            if len - i != 1
               print(io, "^")
               print(io, len - i)
            end
         end
      end
      c = x.data.coeffs[1]
      bracket = isa(c, Poly) && c.data.length != 1
      if c != 0
         if len != 1
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

function -{T <: Ring, S}(a::Poly{T, S})
   len = a.data.length
   z = Poly{T, S}(Array(T, len))
   for i = 1:len
      z.data.coeffs[i] = -a.data.coeffs[i]
   end
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

function +{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   lena = a.data.length
   lenb = b.data.length
   lenz = max(lena, lenb)
   z = Poly{T, S}(Array(T, lenz))
   i = 1

   while i <= min(lena, lenb)
      z.data.coeffs[i] = a.data.coeffs[i] + b.data.coeffs[i]
      i += 1
   end

   while i <= lena
      z.data.coeffs[i] = a.data.coeffs[i]
      i += 1
   end

   while i <= lenb
      z.data.coeffs[i] = b.data.coeffs[i]
      i += 1
   end

   z.data.length = normalise(z, i - 1)

   return z
end
  
function -{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   lena = a.data.length
   lenb = b.data.length
   lenz = max(lena, lenb)
   z = Poly{T, S}(Array(T, lenz))
   i = 1

   while i <= min(lena, lenb)
      z.data.coeffs[i] = a.data.coeffs[i] - b.data.coeffs[i]
      i += 1
   end

   while i <= lena
      z.data.coeffs[i] = a.data.coeffs[i]
      i += 1
   end

   while i <= lenb
      z.data.coeffs[i] = -b.data.coeffs[i]
      i += 1
   end

   z.data.length = normalise(z, i - 1)

   return z
end

function *{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
   lena = a.data.length
   lenb = b.data.length

   if lena == 0 || lenb == 0
      return Poly{T, S}()
   end

   lenz = lena + lenb - 1
   z = Poly{T, S}(Array(T, lenz))

   for i = 1:lena
      z.data.coeffs[i] = a.data.coeffs[i]*b.data.coeffs[1]
   end

   for i = 2:lenb
      z.data.coeffs[lena + i - 1] = a.data.coeffs[lena]*b.data.coeffs[i]
   end

   for i = 1:lena - 1
      for j = 2:lenb
         z.data.coeffs[i + j - 1] += a.data.coeffs[i]*b.data.coeffs[j]
      end
   end
        
   z.data.length = normalise(z, lenz)

   return z
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

function *{T <: Ring, S}(a::Int, b::Poly{T, S})
   len = b.data.length
   z = Poly{T, S}(Array(T, len))
   for i = 1:len
      z.data.coeffs[i] = a*b.data.coeffs[i]
   end
   z.data.length = len
   return z
end

function *{T <: Ring, S}(a::ZZ, b::Poly{T, S})
   len = b.data.length
   z = Poly{T, S}(Array(T, len))
   for i = 1:len
      z.data.coeffs[i] = a*b.data.coeffs[i]
   end
   z.data.length = len
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

function -{S}(x::Int, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_si_sub, :libflint), Void, 
                (Ptr{fmpz_poly}, Int, Ptr{fmpz_poly}), 
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

function truncate{T <: Ring, S}(a::Poly{T, S}, n::Int)
   n < 0 && throw(DomainError())
   
   lena = a.data.length

   if lena <= n
      return a
   end

   lenz = min(lena, n)
   z = Poly{T, S}(Array(T, lenz))

   for i = 1:lenz
      z.data.coeffs[i] = a.data.coeffs[i]
   end

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

function mullow{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S}, n::Int)
   return truncate(a * b, n)
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

function ^{T <: Ring, S}(a::Poly{T, S}, b::Int)
   b < 0 && throw(DomainError())
   # special case powers of x for constructing polynomials efficiently
   if a.data.length == 2 && a.data.coeffs[1] == 0 && a.data.coeffs[2] == 1
      z = Poly{T, S}(Array(T, b + 1))
      z.data.coeffs[b + 1] = a.data.coeffs[2]
      for i = 1:b
         z.data.coeffs[i] = a.data.coeffs[1]
      end
      z.data.length = b + 1
      return z
   elseif a.data.length == 0
      return zero(Poly{T, S})
   elseif a.data.length == 1
      return Poly{T, S}([a.data.coeffs[1]^b])
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
#   Comparisons
#
###########################################################################################

function =={S}(x::Poly{ZZ, S}, y::Int) 
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

=={T <: Ring, S}(x::Poly{T, S}, y::Int) = ((x.data.length == 0 && y == 0)
                                        || (x.data.length == 1 && x.data.coeffs[1] == y))

=={T <: Ring, S}(x::Poly{T, S}, y::ZZ) = ((x.data.length == 0 && y == 0)
                                        || (x.data.length == 1 && x.data.coeffs[1] == y))

=={S}(x::Poly{ZZ, S}, y::Poly{ZZ, S}) = ccall((:fmpz_poly_equal, :libflint), Bool, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}), &(x.data), &(y.data))

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

###########################################################################################
#
#   Division
#
###########################################################################################

function divexact{T <: Ring, S}(a::Poly{T, S}, b::T)
   b == 0 && throw(DivideError())
   z = Poly{T, S}(Array(T, a.data.length))
   for i = 1:a.data.length
      z.data.coeffs[i] = divexact(a.data.coeffs[i], b)
   end
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

function divexact{T <: Ring, S}(f::Poly{T, S}, g::Poly{T, S})
   g == 0 && throw(DivideError())
   if f == 0
      return zero(Poly{T, S})
   end
   lenq = f.data.length - g.data.length + 1
   q = Poly{T, S}(Array(T, lenq))
   for i = 1:lenq
      q.data.coeffs[i] = zero(T)
   end
   x = gen(Poly{T, S})
   leng = g.data.length
   while f.data.length >= leng
      lenf = f.data.length
      q.data.coeffs[lenf - leng + 1] = divexact(f.data.coeffs[lenf], g.data.coeffs[leng])
      f = f - q.data.coeffs[lenf - leng + 1]*g*x^(lenf - leng)
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
   b = g.data.coeffs[g.data.length]
   x = gen(Poly{T, S})
   while f.data.length >= g.data.length
      f = f*b - f.data.coeffs[f.data.length]*g*x^(f.data.length - g.data.length)
   end
   return f
end

function pseudodivrem{T <: Ring, S}(f::Poly{T, S}, g::Poly{T, S})
   g == 0 && throw(DivideError())
   if f.data.length < g.data.length
      return zero(Poly{T, S}), f
   end
   lenq = f.data.length - g.data.length + 1
   q = Poly{T, S}(Array(T, lenq))
   for i = 1:lenq
      q.data.coeffs[i] = zero(T)
   end
   b = g.data.coeffs[g.data.length]
   x = gen(Poly{T, S})
   while f.data.length >= g.data.length
      for i = f.data.length - g.data.length + 2:lenq
         q.data.coeffs[i] *= b
      end
      q.data.coeffs[f.data.length - g.data.length + 1] = f.data.coeffs[f.data.length]
      f = f*b - f.data.coeffs[f.data.length]*g*x^(f.data.length - g.data.length)
   end
   q.data.length = normalise(q, lenq)
   return q, f
end

###########################################################################################
#
#   Content, primitive part and GCD
#
###########################################################################################

function gcd{T <: Ring, S}(a::Poly{T, S}, b::Poly{T, S})
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
      (a, b) = (pseudorem(b, a), a)
   end
   return g*primpart(b)
end

function gcd{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_gcd, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function content{T <: Ring, S}(a::Poly{T, S})
   z = a.data.coeffs[1]
   for i = 2:a.data.length
      z = gcd(z, a.data.coeffs[i])
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

function compose{S}(x::Poly{ZZ, S}, y::Poly{ZZ, S})
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_compose, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

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

function deriv{T <: Ring, S}(a::Poly{T, S})
   if a == 0
      return zero(Poly{T, S})
   end
   len = a.data.length
   z = Poly{T, S}(Array(T, len - 1))
   for i = 1:len - 1
      z.data.coeffs[i] = i*a.data.coeffs[i + 1]
   end
   z.data.length = normalise(z, len - 1)
   return z
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

   Base.convert(::Type{T1}, x::T) = T1([x])
   Base.promote_rule(::Type{T1}, ::Type{T}) = T1

   P = T2.parameters
   while length(P) > 0
      T2 = P[1]
      Base.convert(::Type{T1}, x::T2) = T1([convert(T, x)])
      Base.promote_rule(::Type{T1}, ::Type{T2}) = T1
      P = T2.parameters
   end

   Base.convert(::Type{T1}, x::Integer) = T1([convert(T, x)])
   Base.promote_rule{R <: Integer}(::Type{T1}, ::Type{R}) = T1

   # (Type, gen) 

   return (Poly{T, S}, Poly{T, S}([T(0), T(1)]))
end
