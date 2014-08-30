###########################################################################################
#
#   PowerSeries2.jl : Power series over fields
#
###########################################################################################    

import Rings: PowerSeries, PowerSeriesRing, max, min, isless, Precision, valuation, O,
              coeff, truncate, divexact

export PowerSeries, coeff, truncate, divexact

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

function PowerSeries{S}(::Type{PowerSeries{QQ, S}}, a :: Array{QQ, 1}, n::Precision)
   z = PowerSeries{QQ, S}(PowerSeries, n)
   ccall((:fmpq_poly_init2, :libflint), Void, (Ptr{PowerSeries}, Int), &z, length(a))
   for i = 1:length(a)
      ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, (Ptr{PowerSeries}, Int, Ptr{Fraction}),
            &z, i - 1, &a[i])
   end
   return z
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################    
   
length{S}(x::PowerSeries{QQ, S}) = ccall((:fmpq_poly_length, :libflint), Int, (Ptr{PowerSeries},), &x)

function coeff{S}(x::PowerSeries{QQ, S}, n::Int)
   z = QQ()
   ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, (Ptr{Fraction}, Ptr{PowerSeries}, Int), &z, &x, n)
   return z
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{S}(io::IO, x::PowerSeries{QQ, S})
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fmpq_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{PowerSeries}, Ptr{Uint8}), &x, bytestring(string(S)))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
   if x.prec != nothing
      print(io, "+O(", string(S), "^", x.prec, ")")
   end
end

###########################################################################################
#
#   Unary operators
#
###########################################################################################

function -{S}(x::PowerSeries{QQ, S})
   z = PowerSeries{QQ, S}()
   ccall((:fmpq_poly_neg, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}), 
               &z, &x)
   z.prec = x.prec
   return z
end

###########################################################################################
#
#   Binary operators
#
###########################################################################################

function +{S}(a::PowerSeries{QQ, S}, b::PowerSeries{QQ, S})
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = PowerSeries{QQ, S}()
   z.prec = prec
   ccall((:fmpq_poly_add_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
   return z
end

function -{S}(a::PowerSeries{QQ, S}, b::PowerSeries{QQ, S})
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = PowerSeries{QQ, S}()
   z.prec = prec
   ccall((:fmpq_poly_sub_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
   return z
end

function *{S}(a::PowerSeries{QQ, S}, b::PowerSeries{QQ, S})
   lena = length(a)
   lenb = length(b)
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      return PowerSeries(PowerSeries{QQ, S}, Array(QQ, 0), prec)
   end

   lenz = prec == nothing ? lena + lenb - 1 : min(lena + lenb - 1, prec)

   z = PowerSeries{QQ, S}()
   z.prec = prec
   ccall((:fmpq_poly_mullow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
   return z
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function setcoeff!{S}(z::PowerSeries{QQ, S}, n::Int, x::QQ)
   ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, 
                (Ptr{PowerSeries}, Int, Ptr{Fraction}), 
               &z, n, &x)
end

function mul!{S}(z::PowerSeries{QQ, S}, a::PowerSeries{QQ, S}, b::PowerSeries{QQ, S})
   lena = length(a)
   lenb = length(b)
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   lenz = prec == nothing ? lena + lenb - 1 : min(lena + lenb - 1, prec)
   if lenz < 0
      lenz = 0
   end

   z.prec = prec
   ccall((:fmpq_poly_mullow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
end

function addeq!{S}(a::PowerSeries{QQ, S}, b::PowerSeries{QQ, S},)
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   a.prec = prec
   ccall((:fmpq_poly_add_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &a, &a, &b, lenz)
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{S}(x::Int, y::PowerSeries{QQ, S})
   z = PowerSeries{QQ, S}()
   z.prec = y.prec
   ccall((:fmpq_poly_scalar_mul_si, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &y, x)
   return z
end

function *{S}(x::ZZ, y::PowerSeries{QQ, S})
   z = PowerSeries{QQ, S}()
   z.prec = y.prec
   ccall((:fmpq_poly_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{ZZ}), 
               &z, &y, &x)
   return z
end

function *{S}(x::QQ, y::PowerSeries{QQ, S})
   z = PowerSeries{QQ, S}()
   z.prec = y.prec
   ccall((:fmpq_poly_scalar_mul_fmpq, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{Fraction}), 
               &z, &y, &x)
   return z
end

###########################################################################################
#
#   Comparison
#
###########################################################################################

function =={S}(x::PowerSeries{QQ, S}, y::PowerSeries{QQ, S})
   prec = min(x.prec, y.prec)
   
   n = max(length(x), length(y))
   n = min(n, prec)
   
   return bool(ccall((:fmpq_poly_equal_trunc, :libflint), Cint, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &x, &y, n))
end

###########################################################################################
#
#   Shifting
#
###########################################################################################

function shift_left{S}(x::PowerSeries{QQ, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   z = PowerSeries{QQ, S}()
   z.prec = x.prec + len
   ccall((:fmpq_poly_shift_left, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, len)
   return z
end

function shift_right{S}(x::PowerSeries{QQ, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   if len >= xlen
      return PowerSeries(PowerSeries{QQ, S}, Array(QQ, 0), max(0, x.prec - len))
   end
   z = PowerSeries{QQ, S}()
   z.prec = x.prec - len
   ccall((:fmpq_poly_shift_right, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, len)
   return z
end

###########################################################################################
#
#   Truncation
#
###########################################################################################

function truncate{S}(x::PowerSeries{QQ, S}, prec::Precision)
   prec < 0 && throw(DomainError())
   if x.prec <= prec
      return x
   end
   z = PowerSeries{QQ, S}()
   z.prec = prec
   ccall((:fmpq_poly_set_trunc, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, prec)
   return z
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S}(x::PowerSeries{QQ, S}, y::PowerSeries{QQ, S})
   y == 0 && throw(DivideError())
   v2 = valuation(y)
   v1 = valuation(x)
   if v2 != 0
      if v1 >= v2
         x = shift_right(x, v2)
         y = shift_right(y, v2)
      end
   end
   !isunit(y) && error("Unable to invert power series")
   prec = min(x.prec, y.prec - v2 + v1)
   z = PowerSeries{QQ, S}()
   z.prec = prec
   ccall((:fmpq_poly_div_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, &y, prec)
   return z
end

function divexact{S}(x::PowerSeries{QQ, S}, y::Int)
   y == 0 && throw(DivideError())
   z = PowerSeries{QQ, S}()
   z.prec = x.prec
   ccall((:fmpq_poly_scalar_div_si, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, y)
   return z
end

function divexact{S}(x::PowerSeries{QQ, S}, y::ZZ)
   y == 0 && throw(DivideError())
   z = PowerSeries{QQ, S}()
   z.prec = x.prec
   ccall((:fmpq_poly_scalar_div_fmpz, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{ZZ}), 
               &z, &x, &y)
   return z
end

function divexact{S}(x::PowerSeries{QQ, S}, y::QQ)
   y == 0 && throw(DivideError())
   z = PowerSeries{QQ, S}()
   z.prec = x.prec
   ccall((:fmpq_poly_scalar_div_fmpz, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{Fraction}), 
               &z, &x, &y)
   return z
end

/{S}(x::PowerSeries{QQ, S}, y::Int) = divexact(x, y)

/{S}(x::PowerSeries{QQ, S}, y::ZZ) = divexact(x, y)

/{S}(x::PowerSeries{QQ, S}, y::QQ) = divexact(x, y)

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{S}(a::PowerSeries{QQ, S})
   a == 0 && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   if a.prec == nothing
      a1 = coeff(a, 0)
      length(a) != 1 && error("Unable to invert infinite precision power series")
      return PowerSeries(PowerSeries{QQ, S}, [divexact(QQ(1), a1)], nothing)
   end
   ainv = PowerSeries(PowerSeries{QQ, S}, Array(QQ, 0), a.prec)
   ccall((:fmpq_poly_inv_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &ainv, &a, a.prec)
   return ainv
end

###########################################################################################
#
#   Special functions
#
###########################################################################################

function exp{S}(a::PowerSeries{QQ, S})
   if a == 0
      return PowerSeries(PowerSeries{QQ, S}, [QQ(1)], nothing)
   elseif a.prec == nothing
      error("Unable to compute exponential of infinite precision power series")
   end
   b = PowerSeries(PowerSeries{QQ, S}, Array(QQ, 0), a.prec)
   ccall((:fmpq_poly_exp_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &b, &a, a.prec)
   return b
end