###########################################################################################
#
#   PowerSeries2.jl : Power series over fields
#
###########################################################################################    

import Rings: PowerSeries, PowerSeriesRing, max, min, isless, Precision, valuation, O,
              initps, coeff, truncate

export PowerSeries, coeff, truncate

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

function PowerSeries{S}(::Type{PowerSeries{QQ, S}}, a :: Array{QQ, 1}, n::Precision)
   z = PowerSeries{QQ, S}(initps(), n)
   ccall((:fmpq_poly_init2, :libflint), Void, (Ptr{PowerSeries}, Int), &z, length(a))
   for i = 1:length(a)
      ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, (Ptr{PowerSeries}, Int, Ptr{fmpq}),
            &z, i - 1, &(a[i].data))
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
   ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, (Ptr{fmpq}, Ptr{PowerSeries}, Int), &(z.data), &x, n)
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
   ccall((:fmpq_poly_add, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}), 
               &z, &a, &b)
   ccall((:fmpq_poly_truncate, :libflint), Void, 
                (Ptr{PowerSeries}, Int), 
               &z, lenz)
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
   ccall((:fmpq_poly_sub, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}), 
               &z, &a, &b)
   ccall((:fmpq_poly_truncate, :libflint), Void, 
                (Ptr{PowerSeries}, Int), 
               &z, lenz)
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
      return PowerSeries(PowerSeries{T, S}, Array(T, 0), prec)
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
                (Ptr{PowerSeries}, Int, Ptr{fmpq}), 
               &z, n, &(x.data))
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
   ccall((:fmpq_poly_add, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}), 
               &a, &a, &b)
   ccall((:fmpq_poly_truncate, :libflint), Void, 
                (Ptr{PowerSeries}, Int), 
               &a, lenz)
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
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{fmpq}), 
               &z, &y, &(x.data))
   return z
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
      return PowerSeries(PowerSeries{T, S}, Array(QQ, 0), max(0, x.prec - len))
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
