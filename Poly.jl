export Poly, PolynomialRing, coeff, zero, one, gen, is_zero, is_one, is_gen

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

type Poly{T <: Ring, Symbol} <: Ring
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
         d = new(PolyStruct(a, length(a)))
      end
      return d
   end   

   Poly() = Poly{T, Symbol}(Array(T, 0))
   Poly(a::Integer) = a == 0 ? Poly{T, Symbol}(Array(T, 0)) : Poly{T, Symbol}([T(a)])
   Poly(a::T) = Poly{T, Symbol}([a])
   Poly(a::Poly{T, Symbol}) = a
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

is_zero{S}(x::Poly{ZZ, S}) = x.data.length == 0

is_zero{T <: Ring, S}(a::Poly{T, S}) = a.data.length == 0

is_one{S}(x::Poly{ZZ, S}) = ccall((:__fmpz_poly_is_one, :libflint), Int, (Ptr{fmpz_poly},), &(x.data))

is_one{T <: Ring, S}(a::Poly{T, S}) = a.data.length == 1 && a.data.coeffs[1] == 1

is_gen{S}(x::Poly{ZZ, S}) = ccall((:__fmpz_poly_is_x, :libflint), Int, (Ptr{fmpz_poly},), &(x.data))

is_gen{T <: Ring, S}(a::Poly{T, S}) = a.data.length == 2 && a.data.coeffs[1] == 0 && a.data.coeffs[2] == 1

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
      print(io, "0")
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
#   Powering
#
###########################################################################################


function ^{S}(x::Poly{ZZ, S}, y::Int)
   z = Poly{ZZ, S}()
   ccall((:fmpz_poly_pow, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

function ^{T <: Ring, S}(a::Poly{T, S}, b::Int)
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
      return Poly{T, S}()
   elseif a.data.length == 1
      return Poly{T, S}([a.data.coeffs[0]^b])
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

=={T <: Ring, S}(x::Poly{T, S}, y::Int) = ((x.data.length == 0 && y == 0)
                                        || (x.data.length == 1 && x.data.coeffs[1] == y))

###########################################################################################
#
#   PolynomialRing constructor
#
###########################################################################################

function PolynomialRing(T, s)
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
