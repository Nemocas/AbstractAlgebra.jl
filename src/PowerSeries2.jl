###########################################################################################
#
#   PowerSeries2.jl : Power series over fields
#
###########################################################################################    

import Nemo.Rings: PowerSeries, PowerSeriesRing, max, min, isless, Precision, valuation, O,
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
   finalizer(z, _fmpq_poly_clear_fn)
   return z
end

function PowerSeries{S, T}(::Type{PowerSeries{FinFieldElem{T}, S}}, a :: Array{FinFieldElem{T}, 1}, n::Precision)
   z = PowerSeries{FinFieldElem{T}, S}(PowerSeries, n)
   ctx = FinFieldCtx[T]::fq_ctx
   ccall((:fq_poly_init2, :libflint), Void, (Ptr{PowerSeries}, Int, Ptr{fq_ctx}), &z, length(a), &ctx)
   for i = 1:length(a)
      ccall((:fq_poly_set_coeff, :libflint), Void, (Ptr{PowerSeries}, Int, Ptr{FinFieldElem}, Ptr{fq_ctx}),
            &z, i - 1, &a[i], &ctx)
   end
   finalizer(z, _fq_poly_clear_fn)
   return z
end

function _fmpq_poly_clear_fn{S}(a :: PowerSeries{QQ, S})
   ccall((:fmpq_poly_clear, :libflint), Void, (Ptr{PowerSeries},), &a)
end

function _fq_poly_clear_fn{T, S}(a :: PowerSeries{FinFieldElem{T}, S})
   ccall((:fq_poly_clear, :libflint), Void, (Ptr{PowerSeries},), &a) # clear function doesn't use ctx
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################    
   
length{S}(x::PowerSeries{QQ, S}) = ccall((:fmpq_poly_length, :libflint), Int, (Ptr{PowerSeries},), &x)

length{S, T}(x::PowerSeries{FinFieldElem{T}, S}) = ccall((:fq_poly_length, :libflint), Int, (Ptr{PowerSeries}, Ptr{fq_ctx}), &x, &(FinFieldCtx[T]::fq_ctx))

function coeff{S}(x::PowerSeries{QQ, S}, n::Int)
   if n < 0
      return QQ(0)
   end
   z = QQ()
   ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, (Ptr{Fraction}, Ptr{PowerSeries}, Int), &z, &x, n)
   return z
end

function coeff{S, T}(x::PowerSeries{FinFieldElem{T}, S}, n::Int)
   if n < 0
      return FinFieldElem{T}(0)
   end
   z = FinFieldElem{T}()
   ccall((:fq_poly_get_coeff, :libflint), Void, (Ptr{FinFieldElem}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), &z, &x, n, &(FinFieldCtx[T]::fq_ctx))
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

function show{S, T}(io::IO, x::PowerSeries{FinFieldElem{T}, S})
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fq_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{PowerSeries}, Ptr{Uint8}, Ptr{fq_ctx}), &x, bytestring(string(S)), &(FinFieldCtx[T]::fq_ctx))

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

function -{S, T}(x::PowerSeries{FinFieldElem{T}, S})
   z = PowerSeries{FinFieldElem{T}, S}()
   ccall((:fq_poly_neg, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{fq_ctx}), 
               &z, &x, &(FinFieldCtx[T]::fq_ctx))
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

function +{S, T}(a::PowerSeries{FinFieldElem{T}, S}, b::PowerSeries{FinFieldElem{T}, S})
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = PowerSeries{FinFieldElem{T}, S}()
   z.prec = prec
   ccall((:fq_poly_add_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &z, &a, &b, lenz, &(FinFieldCtx[T]::fq_ctx))
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

function -{S, T}(a::PowerSeries{FinFieldElem{T}, S}, b::PowerSeries{FinFieldElem{T}, S})
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = PowerSeries{FinFieldElem{T}, S}()
   z.prec = prec
   ccall((:fq_poly_sub_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &z, &a, &b, lenz, &(FinFieldCtx[T]::fq_ctx))
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

function *{S, T}(a::PowerSeries{FinFieldElem{T}, S}, b::PowerSeries{FinFieldElem{T}, S})
   lena = length(a)
   lenb = length(b)
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      return PowerSeries(PowerSeries{FinFieldElem{T}, S}, Array(FinFieldElem{T}, 0), prec)
   end

   lenz = prec == nothing ? lena + lenb - 1 : min(lena + lenb - 1, prec)

   z = PowerSeries{FinFieldElem{T}, S}()
   z.prec = prec
   ccall((:fq_poly_mullow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &z, &a, &b, lenz, &(FinFieldCtx[T]::fq_ctx))
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

function setcoeff!{S, T}(z::PowerSeries{FinFieldElem{T}, S}, n::Int, x::FinFieldElem{T})
   ccall((:fq_poly_set_coeff, :libflint), Void, 
                (Ptr{PowerSeries}, Int, Ptr{FinFieldElem{T}}, Ptr{fq_ctx}), 
               &z, n, &x, &(FinFieldCtx[T]::fq_ctx))
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

function mul!{S, T}(z::PowerSeries{FinFieldElem{T}, S}, a::PowerSeries{FinFieldElem{T}, S}, b::PowerSeries{FinFieldElem{T}, S})
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
   ccall((:fq_poly_mullow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &z, &a, &b, lenz, &(FinFieldCtx[T]::fq_ctx))
end

function addeq!{S}(a::PowerSeries{QQ, S}, b::PowerSeries{QQ, S})
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

function addeq!{S, T}(a::PowerSeries{FinFieldElem{T}, S}, b::PowerSeries{FinFieldElem{T}, S})
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   a.prec = prec
   ccall((:fq_poly_add_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &a, &a, &b, lenz, &(FinFieldCtx[T]::fq_ctx))
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

function *{S, T}(x::FinFieldElem{T}, y::PowerSeries{FinFieldElem{T}, S})
   z = PowerSeries{FinFieldElem{T}, S}()
   z.prec = y.prec
   ccall((:fq_poly_scalar_mul_fq, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{FinFieldElem{T}}, Ptr{fq_ctx}), 
               &z, &y, &x, &(FinFieldCtx[T]::fq_ctx))
   return z
end

*{S, T}(x::Int, y::PowerSeries{FinFieldElem{T}, S}) = FinFieldElem{T}(x)*y

*{S, T}(x::ZZ, y::PowerSeries{FinFieldElem{T}, S}) = FinFieldElem{T}(x)*y

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

function =={S, T}(x::PowerSeries{FinFieldElem{T}, S}, y::PowerSeries{FinFieldElem{T}, S})
   prec = min(x.prec, y.prec)
   
   n = max(length(x), length(y))
   n = min(n, prec)
   
   return bool(ccall((:fq_poly_equal_trunc, :libflint), Cint, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &x, &y, n, &(FinFieldCtx[T]::fq_ctx)))
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

function shift_left{S, T}(x::PowerSeries{FinFieldElem{T}, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   z = PowerSeries{FinFieldElem{T}, S}()
   z.prec = x.prec + len
   ccall((:fq_poly_shift_left, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &z, &x, len, &(FinFieldCtx[T]::fq_ctx))
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

function shift_right{S, T}(x::PowerSeries{FinFieldElem{T}, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   if len >= xlen
      return PowerSeries(PowerSeries{FinFieldElem{T}, S}, Array(FinFieldElem{T}, 0), max(0, x.prec - len))
   end
   z = PowerSeries{FinFieldElem{T}, S}()
   z.prec = x.prec - len
   ccall((:fq_poly_shift_right, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &z, &x, len, &(FinFieldCtx[T]::fq_ctx))
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

function truncate{S, T}(x::PowerSeries{FinFieldElem{T}, S}, prec::Precision)
   prec < 0 && throw(DomainError())
   if x.prec <= prec
      return x
   end
   z = PowerSeries{FinFieldElem{T}, S}()
   z.prec = prec
   ccall((:fq_poly_set_trunc, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &z, &x, prec, &(FinFieldCtx[T]::fq_ctx))
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
   if prec == nothing
      length(y) != 1 && error("Unable to invert infinite precision power series")
      return divexact(x, coeff(y, 0))
   end
   z = PowerSeries{QQ, S}()
   z.prec = prec
   ccall((:fmpq_poly_div_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, &y, prec)
   return z
end

function divexact{S, T}(x::PowerSeries{FinFieldElem{T}, S}, y::PowerSeries{FinFieldElem{T}, S})
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
   if prec == nothing
      length(y) != 1 && error("Unable to invert infinite precision power series")
      return divexact(x, coeff(y, 0))
   end
   z = PowerSeries{FinFieldElem{T}, S}()
   z.prec = prec
   ccall((:fq_poly_div_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &z, &x, &y, prec, &(FinFieldCtx[T]::fq_ctx))
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
   ccall((:fmpq_poly_scalar_div_fmpq, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{Fraction}), 
               &z, &x, &y)
   return z
end

function divexact{S, T}(x::PowerSeries{FinFieldElem{T}, S}, y::FinFieldElem{T})
   y == 0 && throw(DivideError())
   z = PowerSeries{FinFieldElem{T}, S}()
   z.prec = x.prec
   ccall((:fq_poly_scalar_div_fq, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{FinFieldElem}, Ptr{fq_ctx}), 
               &z, &x, &y, &(FinFieldCtx[T]::fq_ctx))
   return z
end

divexact{S, T}(x::PowerSeries{FinFieldElem{T}, S}, y::ZZ) = divexact(x, FinFieldElem{T}(y))

divexact{S, T}(x::PowerSeries{FinFieldElem{T}, S}, y::Int) = divexact(x, FinFieldElem{T}(y))

/{S}(x::PowerSeries{QQ, S}, y::Int) = divexact(x, y)

/{S}(x::PowerSeries{QQ, S}, y::ZZ) = divexact(x, y)

/{S}(x::PowerSeries{QQ, S}, y::QQ) = divexact(x, y)

/{S, T}(x::PowerSeries{FinFieldElem{T}, S}, y::Int) = divexact(x, y)

/{S, T}(x::PowerSeries{FinFieldElem{T}, S}, y::ZZ) = divexact(x, y)

/{S, T}(x::PowerSeries{FinFieldElem{T}, S}, y::FinFieldElem{T}) = divexact(x, y)

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

function inv{S, T}(a::PowerSeries{FinFieldElem{T}, S})
   a == 0 && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   if a.prec == nothing
      a1 = coeff(a, 0)
      length(a) != 1 && error("Unable to invert infinite precision power series")
      return PowerSeries(PowerSeries{FinFieldElem{T}, S}, [divexact(FinFieldElem{T}(1), a1)], nothing)
   end
   ainv = PowerSeries(PowerSeries{FinFieldElem{T}, S}, Array(FinFieldElem{T}, 0), a.prec)
   ccall((:fq_poly_inv_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Ptr{fq_ctx}), 
               &ainv, &a, a.prec, &(FinFieldCtx[T]::fq_ctx))
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