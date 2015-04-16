###########################################################################################
#
#   PowerSeries.jl : Power series over rings
#
###########################################################################################    

export PowerSeries, PowerSeriesRing, O, valuation, min, max, isless, Precision, exp

import Base: min, max, isless, exp

###########################################################################################
#
#   Precision type
#
###########################################################################################

typealias Precision Union(Nothing, Int)

max(a::Nothing, b::Int) = nothing

max(a::Int, b::Nothing) = nothing

max(a::Nothing, b::Nothing) = nothing

min(a::Nothing, b::Int) = b

min(a::Int, b::Nothing) = a

min(a::Nothing, b::Nothing) = nothing

isless(a::Nothing, b::Int) = false

isless(a::Int, b::Nothing) = true

isless(a::Nothing, b::Nothing) = false

+(a::Nothing, b::Int) = nothing

+(a::Int, b::Nothing) = nothing

+(a::Nothing, b::Nothing) = nothing

-(a::Nothing, b::Int) = nothing

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

type PowerSeries{T <: Ring, S} <: Ring
   arr::Ptr{Void}
   length::Int
   alloc::Int
   inv::Int
   prec :: Precision
   coeffs :: Array{T, 1}
   
   PowerSeries(a :: Array{T, 1}, n :: Precision) = new(C_NULL, 0, 0, 0, n, a)   

   PowerSeries(::Type{PowerSeries}, n :: Precision) = new(C_NULL, 0, 0, 0, n)
   
   PowerSeries() = PowerSeries(PowerSeries{T, S}, Array(T, 0), nothing)
   
   PowerSeries(a::Integer) = a == 0 ? PowerSeries(PowerSeries{T, S}, Array(T, 0), nothing) : PowerSeries(PowerSeries{T, S}, [T(a)], nothing)
   PowerSeries(a::T) = PowerSeries(PowerSeries{T, S}, [a], nothing)
   PowerSeries(a::PowerSeries{T, S}) = a
   PowerSeries{R <: Ring}(a::R) = convert(PowerSeries{T, S}, a)
end

function PowerSeries{T, S}(::Type{PowerSeries{T, S}}, a :: Array{T, 1}, n :: Precision)
   z = PowerSeries{T, S}(a, n)
   z.length = normalise(z, length(a))
   return z
end

function PowerSeries{S}(::Type{PowerSeries{ZZ, S}}, a :: Array{ZZ, 1}, n::Precision)
   z = PowerSeries{ZZ, S}(PowerSeries, n)
   ccall((:fmpz_poly_init2, :libflint), Void, (Ptr{PowerSeries}, Int), &z, length(a))
   for i = 1:length(a)
      ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Void, (Ptr{PowerSeries}, Int, Ptr{ZZ}),
            &z, i - 1, &a[i])
   end
   finalizer(z, _fmpz_poly_clear_fn)
   return z
end

function PowerSeries{S, M}(::Type{PowerSeries{Residue{ZZ, M}, S}}, a :: Array{Residue{ZZ, M}, 1}, n::Precision)
   z = PowerSeries{Residue{ZZ, M}, S}(PowerSeries, n)
   ccall((:fmpz_mod_poly_init2, :libflint), Void, (Ptr{PowerSeries}, Ptr{ZZ}, Int), &z, &(ResidueModulus[M]::ZZ), length(a))
   for i = 1:length(a)
      ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Void, (Ptr{PowerSeries}, Int, Ptr{ZZ}),
            &z, i - 1, &(a[i].data))
   end
   finalizer(z, _fmpz_mod_poly_clear_fn)
   return z
end

function _fmpz_poly_clear_fn{S}(a :: PowerSeries{ZZ, S})
   ccall((:fmpz_poly_clear, :libflint), Void, (Ptr{PowerSeries},), &a)
end
   
function _fmpz_mod_poly_clear_fn{M, S}(a :: PowerSeries{Residue{ZZ, M}, S})
   ccall((:fmpz_mod_poly_clear, :libflint), Void, (Ptr{PowerSeries},), &a)
end


function O{T, S}(a :: PowerSeries{T, S})
   prec = length(a) - 1
   prec < 0 && throw(DivideError())
   a.prec != nothing && error("Invalid power series monomial in O()")
   return PowerSeries(PowerSeries{T, S}, Array(T, 0), prec)
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################    
   
length{T <: Ring, S}(x::PowerSeries{T, S}) = x.length

length{S}(x::PowerSeries{ZZ, S}) = ccall((:fmpz_poly_length, :libflint), Int, (Ptr{PowerSeries},), &x)

length{S, M}(x::PowerSeries{Residue{ZZ, M}, S}) = ccall((:fmpz_mod_poly_length, :libflint), Int, (Ptr{PowerSeries},), &x)

degree{T <: Ring, S}(x::PowerSeries{T, S}) = length(x) - 1

function normalise{T <: Ring, S}(a::PowerSeries{T, S}, len::Int)
   while len > 0 && iszero(a.coeffs[len]) # cannot use coeff(a, len - 1) here
      len -= 1
   end

   return len
end

coeff{T <: Ring, S}(a::PowerSeries{T, S}, n::Int) = n < 0 || n >= a.length ? T(0) : a.coeffs[n + 1]

function coeff{S}(x::PowerSeries{ZZ, S}, n::Int)
   z = ZZ()
   if n < 0
      return 0
   end
   ccall((:fmpz_poly_get_coeff_fmpz, :libflint), Void, (Ptr{ZZ}, Ptr{PowerSeries}, Int), &z, &x, n)
   return z
end

function coeff{S, M}(x::PowerSeries{Residue{ZZ, M}, S}, n::Int)
   z = ZZ()
   if n < 0
      return Residue{ZZ, M}(0)
   end
   ccall((:fmpz_mod_poly_get_coeff_fmpz, :libflint), Void, (Ptr{ZZ}, Ptr{PowerSeries}, Int), &z, &x, n)
   return Residue{ZZ, M}(z)
end

zero{T <: Ring, S}(::Type{PowerSeries{T, S}}) = PowerSeries{T, S}(0)

one{T <: Ring, S}(::Type{PowerSeries{T, S}}) = PowerSeries{T, S}(1)

iszero{T <: Ring, S}(a::PowerSeries{T, S}) = a.prec == nothing && length(a) == 0

isone{T <: Ring, S}(a::PowerSeries{T, S}) = a.prec == nothing && length(a) == 1 && isone(coeff(a, 0))

gen{T <: Ring, S}(::Type{PowerSeries{T, S}}) = PowerSeries(PowerSeries{T, S}, [T(0), T(1)], nothing)

isgen{T <: Ring, S}(a::PowerSeries{T, S}) = a.prec == nothing && length(a) == 2 && iszero(coeff(a, 0)) && isone(coeff(a, 1))

isunit{T <: Ring, S}(a::PowerSeries{T, S}) = isunit(coeff(a, 0))

function valuation{T <: Ring, S}(a::PowerSeries{T, S})
   if length(a) == 0
      return a.prec
   end
   for i = 1:length(a)
      if coeff(a, i - 1) != 0
         return i - 1
      end
   end
   error("Power series is not normalised")
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{T <: Ring, S}(io::IO, x::PowerSeries{T, S})
   len = x.length

   if len == 0
      print(io, zero(T))
   else
      coeff_printed = false
      c = x.coeffs[1]
      bracket = needs_parentheses(c)
      if !iszero(c)
         if bracket
            print(io, "(")
         end
         show(io, c)
         if bracket
            print(io, ")")
         end
         coeff_printed = true
      end
      for i = 2:len
         c = x.coeffs[i]
         bracket = needs_parentheses(c)
         if !iszero(c)
            if coeff_printed && !is_negative(c)
               print(io, "+")
            end
            if !isone(c) && (c != -1 || show_minus_one(typeof(c)))
               if bracket
                  print(io, "(")
               end
               show(io, c)
               if bracket
                  print(io, ")")
               end
               print(io, "*")
            end
            if c == -1 && !show_minus_one(typeof(c))
               print(io, "-")
            end
            print(io, string(S))
            if i != 2
               print(io, "^")
               print(io, i - 1)
            end
            coeff_printed = true
         end
      end
      
   end
   if x.prec != nothing
      print(io, "+O(", string(S), "^", x.prec, ")")
   end
end

function show{S}(io::IO, x::PowerSeries{ZZ, S})
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fmpz_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{PowerSeries}, Ptr{Uint8}), &x, bytestring(string(S)))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
   if x.prec != nothing
      print(io, "+O(", string(S), "^", x.prec, ")")
   end
end

function show{S, M}(io::IO, x::PowerSeries{Residue{ZZ, M}, S})
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fmpz_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{PowerSeries}, Ptr{Uint8}), &x, bytestring(string(S)))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
   if x.prec != nothing
      print(io, "+O(", string(S), "^", x.prec, ")")
   end
end

function show{T <: Ring, S}(io::IO, ::Type{PowerSeries{T, S}})
   print(io, "Univariate power series ring in ", string(S), " over ")
   show(io, T)
end

needs_parentheses{T <: Ring, S}(x::PowerSeries{T, S}) = length(x) > 1

is_negative{T <: Ring, S}(x::PowerSeries{T, S}) = length(x) <= 1 && is_negative(coeff(x, 0))

show_minus_one{T <: Ring, S}(::Type{PowerSeries{T, S}}) = show_minus_one(T)

###########################################################################################
#
#   Unary operators
#
###########################################################################################

function -{T <: Ring, S}(a::PowerSeries{T, S})
   len = a.length
   d = Array(T, len)
   for i = 1:len
      d[i] = -a.coeffs[i]
   end
   z = PowerSeries(PowerSeries{T, S}, d, a.prec)
   z.length = len
   return z
end

function -{S}(x::PowerSeries{ZZ, S})
   z = PowerSeries{ZZ, S}()
   ccall((:fmpz_poly_neg, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}), 
               &z, &x)
   z.prec = x.prec
   return z
end

function -{S, M}(x::PowerSeries{Residue{ZZ, M}, S})
   z = PowerSeries{Residue{ZZ, M}, S}()
   ccall((:fmpz_mod_poly_neg, :libflint), Void, 
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

function +{T <: Ring, S}(a::PowerSeries{T, S}, b::PowerSeries{T, S})
   lena = a.length
   lenb = b.length
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = a.coeffs[i] + b.coeffs[i]
      i += 1
   end

   while i <= lena
      d[i] = a.coeffs[i]
      i += 1
   end

   while i <= lenb
      d[i] = b.coeffs[i]
      i += 1
   end

   z = PowerSeries(PowerSeries{T, S}, d, prec)

   z.length = normalise(z, i - 1)

   return z
end
  
function +{S}(a::PowerSeries{ZZ, S}, b::PowerSeries{ZZ, S})
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = PowerSeries{ZZ, S}()
   z.prec = prec
   ccall((:fmpz_poly_add_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
   return z
end

function +{S, M}(a::PowerSeries{Residue{ZZ, M}, S}, b::PowerSeries{Residue{ZZ, M}, S})
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = PowerSeries{Residue{ZZ, M}, S}()
   z.prec = prec
   ccall((:fmpz_mod_poly_add_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
   return z
end

function -{T <: Ring, S}(a::PowerSeries{T, S}, b::PowerSeries{T, S})
   lena = a.length
   lenb = b.length
   
   prec = min(a.prec, b.prec)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = a.coeffs[i] - b.coeffs[i]
      i += 1
   end

   while i <= lena
      d[i] = a.coeffs[i]
      i += 1
   end

   while i <= lenb
      d[i] = -b.coeffs[i]
      i += 1
   end

   z = PowerSeries(PowerSeries{T, S}, d, prec)

   z.length = normalise(z, i - 1)

   return z
end

function -{S}(a::PowerSeries{ZZ, S}, b::PowerSeries{ZZ, S})
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = PowerSeries{ZZ, S}()
   z.prec = prec
   ccall((:fmpz_poly_sub_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
   return z
end

function -{S, M}(a::PowerSeries{Residue{ZZ, M}, S}, b::PowerSeries{Residue{ZZ, M}, S})
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = PowerSeries{Residue{ZZ, M}, S}()
   z.prec = prec
   ccall((:fmpz_mod_poly_sub_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
   return z
end

function *{T <: Ring, S}(a::PowerSeries{T, S}, b::PowerSeries{T, S})
   lena = a.length
   lenb = b.length
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      return PowerSeries(PowerSeries{T, S}, Array(T, 0), prec)
   end

   t = T()

   lenz = prec == nothing ? lena + lenb - 1 : min(lena + lenb - 1, prec)

   d = Array(T, lenz)

   for i = 1:min(lena, lenz)
      d[i] = a.coeffs[i]*b.coeffs[1]
   end

   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = a.coeffs[lena]*b.coeffs[j]
      end
   end

   z = PowerSeries(PowerSeries{T, S}, d, prec)

   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            mul!(t, a.coeffs[i], b.coeffs[j])
            addeq!(z.coeffs[i + j - 1], t)
         end
      end
   end
        
   z.length = normalise(z, lenz)

   return z
end

function *{S}(a::PowerSeries{ZZ, S}, b::PowerSeries{ZZ, S})
   lena = length(a)
   lenb = length(b)
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      return PowerSeries(PowerSeries{ZZ, S}, Array(ZZ, 0), prec)
   end

   lenz = prec == nothing ? lena + lenb - 1 : min(lena + lenb - 1, prec)

   z = PowerSeries{ZZ, S}()
   z.prec = prec
   ccall((:fmpz_poly_mullow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
   return z
end

function *{S, M}(a::PowerSeries{Residue{ZZ, M}, S}, b::PowerSeries{Residue{ZZ, M}, S})
   lena = length(a)
   lenb = length(b)
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      return PowerSeries(PowerSeries{Residue{ZZ, M}, S}, Array(Residue{ZZ, M}, 0), prec)
   end

   lenz = prec == nothing ? lena + lenb - 1 : min(lena + lenb - 1, prec)

   z = PowerSeries{Residue{ZZ, M}, S}()
   z.prec = prec
   ccall((:fmpz_mod_poly_mullow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
   return z
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function fit!{T <: Ring, S}(c::PowerSeries{T, S}, n::Int)
   if c.length < n
      t = c.coeffs
      c.coeffs = Array(T, n)
      for i = 1:c.length
         c.coeffs[i] = t[i]
      end
      for i = c.length + 1:n
         c.coeffs[i] = zero(T)
      end
   end
end

function setcoeff!{T <: Ring, S}(c::PowerSeries{T, S}, n::Int, a::T)
   if (a != 0 && (c.prec == nothing || c.prec > n)) || n + 1 <= c.length
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(c.length, n + 1)
      # don't normalise
   end
end

function setcoeff!{S}(z::PowerSeries{ZZ, S}, n::Int, x::ZZ)
   ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Void, 
                (Ptr{PowerSeries}, Int, Ptr{ZZ}), 
               &z, n, &x)
end

function setcoeff!{S, M}(z::PowerSeries{Residue{ZZ, M}, S}, n::Int, x::Residue{ZZ, M})
   ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Void, 
                (Ptr{PowerSeries}, Int, Ptr{ZZ}), 
               &z, n, &(x.data))
end

function mul!{T <: Ring, S}(c::PowerSeries{T, S}, a::PowerSeries{T, S}, b::PowerSeries{T, S})
   lena = a.length
   lenb = b.length

   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      c.length = 0
   else
      t = T()

      lenc = prec == nothing ? lena + lenb - 1 : min(lena + lenb - 1, prec)
      fit!(c, lenc)

      for i = 1:min(lena, lenc)
         mul!(c.coeffs[i], a.coeffs[i], b.coeffs[1])
      end

      if lenc > lena
         for i = 2:min(lenb, lenc - lena + 1)
            mul!(c.coeffs[lena + i - 1], a.coeffs[lena], b.coeffs[i])
         end
      end

      for i = 1:lena - 1
         if lenc > i
            for j = 2:min(lenb, lenc - i + 1)
               mul!(t, a.coeffs[i], b.coeffs[j])
               addeq!(c.coeffs[i + j - 1], t)
            end
         end
      end
        
      c.length = normalise(c, lenc)
   end
   c.prec = prec
end

function mul!{S}(z::PowerSeries{ZZ, S}, a::PowerSeries{ZZ, S}, b::PowerSeries{ZZ, S})
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
   ccall((:fmpz_poly_mullow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
end

function mul!{S, M}(z::PowerSeries{Residue{ZZ, M}, S}, a::PowerSeries{Residue{ZZ, M}, S}, b::PowerSeries{Residue{ZZ, M}, S})
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
   ccall((:fmpz_mod_poly_mullow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, &b, lenz)
end

function addeq!{T <: Ring, S}(c::PowerSeries{T, S}, a::PowerSeries{T, S})
   lenc = c.length
   lena = a.length
   
   prec = min(a.prec, c.prec)
   
   lena = min(lena, prec)
   lenc = min(lenc, prec)

   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      addeq!(c.coeffs[i], a.coeffs[i])
   end
   c.length = normalise(c, len)
   c.prec = prec
end

function addeq!{S}(a::PowerSeries{ZZ, S}, b::PowerSeries{ZZ, S},)
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   a.prec = prec
   ccall((:fmpz_poly_add_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &a, &a, &b, lenz)
end

function addeq!{S, M}(a::PowerSeries{Residue{ZZ, M}, S}, b::PowerSeries{Residue{ZZ, M}, S},)
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   a.prec = prec
   ccall((:fmpz_mod_poly_add_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &a, &a, &b, lenz)
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{T <: Ring, S}(a::Int, b::PowerSeries{T, S})
   len = b.length
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = PowerSeries(PowerSeries{T, S}, d, b.prec)
   z.length = normalise(z, len)
   return z
end

function *{T <: Ring, S}(a::ZZ, b::PowerSeries{T, S})
   len = b.length
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = PowerSeries(PowerSeries{T, S}, d, b.prec)
   z.length = normalise(z, len)
   return z
end

function *{S}(x::Int, y::PowerSeries{ZZ, S})
   z = PowerSeries{ZZ, S}()
   z.prec = y.prec
   ccall((:fmpz_poly_scalar_mul_si, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &y, x)
   return z
end

function *{S}(x::ZZ, y::PowerSeries{ZZ, S})
   z = PowerSeries{ZZ, S}()
   z.prec = y.prec
   ccall((:fmpz_poly_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{ZZ}), 
               &z, &y, &x)
   return z
end

function *{S, M}(x::ZZ, y::PowerSeries{Residue{ZZ, M}, S})
   z = PowerSeries{Residue{ZZ, M}, S}()
   z.prec = y.prec
   ccall((:fmpz_mod_poly_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{ZZ}), 
               &z, &y, &x)
   return z
end

*{S, M}(x::Int, y::PowerSeries{Residue{ZZ, M}, S}) = ZZ(x)*y

*{T <: Ring, S}(a::PowerSeries{T, S}, b::Int) = b*a

*{T <: Ring, S}(a::PowerSeries{T, S}, b::ZZ) = b*a

###########################################################################################
#
#   Shifting
#
###########################################################################################

function shift_left{T <: Ring, S}(x::PowerSeries{T, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = x.length
   v = Array(T, xlen + len)
   for i = 1:len
      v[i] = zero(T)
   end
   for i = 1:xlen
      v[i + len] = coeff(x, i - 1)
   end
   return PowerSeries(PowerSeries{T, S}, v, x.prec + len)
end

function shift_left{S}(x::PowerSeries{ZZ, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   z = PowerSeries{ZZ, S}()
   z.prec = x.prec + len
   ccall((:fmpz_poly_shift_left, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, len)
   return z
end

function shift_left{S, M}(x::PowerSeries{Residue{ZZ, M}, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   z = PowerSeries{Residue{ZZ, M}, S}()
   z.prec = x.prec + len
   ccall((:fmpz_mod_poly_shift_left, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, len)
   return z
end

function shift_right{T <: Ring, S}(x::PowerSeries{T, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = x.length
   if len >= xlen
      return PowerSeries(PowerSeries{T, S}, Array(T, 0), max(0, x.prec - len))
   end
   v = Array(T, xlen - len)
   for i = 1:xlen - len
      v[i] = coeff(x, i + len - 1)
   end
   return PowerSeries(PowerSeries{T, S}, v, x.prec - len)
end

function shift_right{S}(x::PowerSeries{ZZ, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   if len >= xlen
      return PowerSeries(PowerSeries{ZZ, S}, Array(ZZ, 0), max(0, x.prec - len))
   end
   z = PowerSeries{ZZ, S}()
   z.prec = x.prec - len
   ccall((:fmpz_poly_shift_right, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, len)
   return z
end

function shift_right{S, M}(x::PowerSeries{Residue{ZZ, M}, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   if len >= xlen
      return PowerSeries(PowerSeries{Residue{ZZ, M}, S}, Array(Residue{ZZ, M}, 0), max(0, x.prec - len))
   end
   z = PowerSeries{Residue{ZZ, M}, S}()
   z.prec = x.prec - len
   ccall((:fmpz_mod_poly_shift_right, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, len)
   return z
end

###########################################################################################
#
#   Truncation
#
###########################################################################################

function truncate{T<: Ring, S}(x::PowerSeries{T, S}, prec::Precision)
   prec < 0 && throw(DomainError())
   if x.prec <= prec
      return x
   end
   d = Array(T, prec)
   for i = 1:min(prec, x.length)
      d[i] = coeff(x, i - 1)
   end
   for i = x.length + 1:prec
      d[i] = zero(T)
   end
   r = PowerSeries(PowerSeries{T, S}, d, prec)
   r.length = normalise(r, prec)
   return r
end

function truncate{S}(x::PowerSeries{ZZ, S}, prec::Precision)
   prec < 0 && throw(DomainError())
   if x.prec <= prec
      return x
   end
   z = PowerSeries{ZZ, S}()
   z.prec = prec
   ccall((:fmpz_poly_set_trunc, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, prec)
   return z
end

function truncate{S, M}(x::PowerSeries{Residue{ZZ, M}, S}, prec::Precision)
   prec < 0 && throw(DomainError())
   if x.prec <= prec
      return x
   end
   z = PowerSeries{Residue{ZZ, M}, S}()
   z.prec = prec
   ccall((:fmpz_mod_poly_set_trunc, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, prec)
   return z
end

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{T <: Ring, S}(a::PowerSeries{T, S}, b::Int)
   b < 0 && throw(DomainError())
   # special case powers of x for constructing power series efficiently
   if a.prec == nothing && isgen(a)
      d = Array(T, b + 1)
      d[b + 1] = coeff(a, 1)
      for i = 1:b
         d[i] = coeff(a, 0)
      end
      return PowerSeries(PowerSeries{T, S}, d, nothing)
   elseif length(a) == 0
      return PowerSeries(PowerSeries{T, S}, Array(T, 0), a.prec + (b - 1)*valuation(a))
   elseif length(a) == 1
      return PowerSeries(PowerSeries{T, S}, [a.coeffs[1]^b], a.prec)
   elseif b == 0
      return PowerSeries(PowerSeries{T, S}, [T(1)], nothing)
   else
      bit = ~((~Uint(0)) >> 1)
      while (Uint(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = z*z
         if (Uint(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

function ^{S}(a::PowerSeries{ZZ, S}, b::Int)
   b < 0 && throw(DomainError())
   if length(a) == 0
      return PowerSeries(PowerSeries{ZZ, S}, Array(ZZ, 0), a.prec + (b - 1)*valuation(a))
   elseif length(a) == 1
      return PowerSeries(PowerSeries{ZZ, S}, [coeffs(a, 0)^b], a.prec)
   elseif b == 0
      return PowerSeries(PowerSeries{ZZ, S}, [ZZ(1)], nothing)
   elseif a.prec == nothing
      z = PowerSeries(PowerSeries{ZZ, S}, Array(ZZ, 0), nothing)
      ccall((:fmpz_poly_pow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, b)
      return z
   else
      prec = a.prec + (b - 1)*valuation(a)
      z = PowerSeries(PowerSeries{ZZ, S}, Array(ZZ, 0), prec)
      ccall((:fmpz_poly_pow_trunc, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Int), 
               &z, &a, b, prec)
      return z
   end
end

function ^{S, M}(a::PowerSeries{Residue{ZZ, M}, S}, b::Int)
   b < 0 && throw(DomainError())
   if length(a) == 0
      return PowerSeries(PowerSeries{Residue{ZZ, M}, S}, Array(ZZ, 0), a.prec + (b - 1)*valuation(a))
   elseif length(a) == 1
      return PowerSeries(PowerSeries{Residue{ZZ, M}, S}, [coeffs(a, 0)^b], a.prec)
   elseif b == 0
      return PowerSeries(PowerSeries{Residue{ZZ, M}, S}, [ZZ(1)], nothing)
   elseif a.prec == nothing
      z = PowerSeries(PowerSeries{Residue{ZZ, M}, S}, Array(Residue{ZZ, M}, 0), nothing)
      ccall((:fmpz_mod_poly_pow, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &a, b)
      return z
   else
      prec = a.prec + (b - 1)*valuation(a)
      z = PowerSeries(PowerSeries{Residue{ZZ, M}, S}, Array(Residue{ZZ, M}, 0), prec)
      ccall((:fmpz_mod_poly_pow_trunc, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int, Int), 
               &z, &a, b, prec)
      return z
   end
end

###########################################################################################
#
#   Comparisons
#
###########################################################################################

=={T <: Ring, S}(x::PowerSeries{T, S}, y::Int) = x.prec == 0 || ((length(x) == 0 && y == 0)
                                        || (length(x) == 1 && coeff(x, 0) == y))

=={T <: Ring, S}(x::PowerSeries{T, S}, y::ZZ) = x.prec == 0 || ((length(x) == 0 && y == 0)
                                        || (length(x) == 1 && coeff(x, 0) == y))

function =={T<: Ring, S}(x::PowerSeries{T, S}, y::PowerSeries{T, S})
   prec = min(x.prec, y.prec)
   
   m1 = min(length(x), length(y))
   m2 = max(length(x), length(y))
   
   m1 = min(m1, prec)
   m2 = min(m2, prec)
   if length(x) >= m2
      for i = m1 + 1: m2
         if coeff(x, i - 1) != 0
            return false
          end
      end
   else
      for i = m1 + 1: m2
         if coeff(y, i - 1) != 0
            return false
          end
      end
   end
           
   for i = 1:m1
      if coeff(x, i - 1) != coeff(y, i - 1)
         return false
      end
   end

   return true
end

function =={S}(x::PowerSeries{ZZ, S}, y::PowerSeries{ZZ, S})
   prec = min(x.prec, y.prec)
   
   n = max(length(x), length(y))
   n = min(n, prec)
   
   return bool(ccall((:fmpz_poly_equal_trunc, :libflint), Cint, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &x, &y, n))
end

function =={S, M}(x::PowerSeries{Residue{ZZ, M}, S}, y::PowerSeries{Residue{ZZ, M}, S})
   prec = min(x.prec, y.prec)
   
   n = max(length(x), length(y))
   n = min(n, prec)
   
   return bool(ccall((:fmpz_mod_poly_equal_trunc, :libflint), Cint, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &x, &y, n))
end

=={T<: Ring, S}(x::Int, y::PowerSeries{T, S}) = y == x

=={T<: Ring, S}(x::ZZ, y::PowerSeries{T, S}) = y == x

function isequal{T <: Ring, S}(x::PowerSeries{T, S}, y::PowerSeries{T, S})
   if x.prec != y.prec || length(x) != length(y)
      return false
   end
   for i = 1:length(x)
      if !isequal(coeff(x, i - 1), coeff(y, i - 1))
         return false
      end
   end
   return true
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{T<: Ring, S}(x::PowerSeries{T, S}, y::PowerSeries{T, S})
   y == 0 && throw(DivideError())
   v2 = valuation(y)
   if v2 != 0
      v1 = valuation(x)
      if v1 >= v2
         x = shift_right(x, v2)
         y = shift_right(y, v2)
      end
   end
   y = truncate(y, x.prec)
   return x*inv(y)
end

function divexact{S}(x::PowerSeries{ZZ, S}, y::PowerSeries{ZZ, S})
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
   z = PowerSeries{ZZ, S}()
   z.prec = prec
   ccall((:fmpz_poly_div_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, &y, prec)
   return z
end

function divexact{S, M}(x::PowerSeries{Residue{ZZ, M}, S}, y::PowerSeries{Residue{ZZ, M}, S})
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
   z = PowerSeries{Residue{ZZ, M}, S}()
   z.prec = prec
   ccall((:fmpz_mod_poly_div_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, &y, prec)
   return z
end

function divexact{T<: Ring, S}(x::PowerSeries{T, S}, y::Int)
   y == 0 && throw(DivideError())
   lenx = length(x)
   d = Array(T, lenx)
   for i = 1:lenx
      d[i] = divexact(coeff(x, i - 1), y)
   end
   return PowerSeries(PowerSeries{T, S}, d, x.prec)
end

function divexact{T<: Ring, S}(x::PowerSeries{T, S}, y::ZZ)
   y == 0 && throw(DivideError())
   lenx = length(x)
   d = Array(T, lenx)
   for i = 1:lenx
      d[i] = divexact(coeff(x, i - 1), y)
   end
   return PowerSeries(PowerSeries{T, S}, d, x.prec)
end

function divexact{T<: Ring, S}(x::PowerSeries{T, S}, y::T)
   y == 0 && throw(DivideError())
   lenx = length(x)
   d = Array(T, lenx)
   for i = 1:lenx
      d[i] = divexact(coeff(x, i - 1), y)
   end
   return PowerSeries(PowerSeries{T, S}, d, x.prec)
end

function divexact{S}(x::PowerSeries{ZZ, S}, y::Int)
   y == 0 && throw(DivideError())
   z = PowerSeries{ZZ, S}()
   z.prec = x.prec
   ccall((:fmpz_poly_scalar_divexact_si, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &z, &x, y)
   return z
end

function divexact{S}(x::PowerSeries{ZZ, S}, y::ZZ)
   y == 0 && throw(DivideError())
   z = PowerSeries{ZZ, S}()
   z.prec = x.prec
   ccall((:fmpz_poly_scalar_divexact_fmpz, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{ZZ}), 
               &z, &x, &y)
   return z
end

function divexact{S, M}(x::PowerSeries{Residue{ZZ, M}, S}, y::ZZ)
   y == 0 && throw(DivideError())
   z = PowerSeries{Residue{ZZ, M}, S}()
   z.prec = x.prec
   ccall((:fmpz_mod_poly_scalar_div_fmpz, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Ptr{ZZ}), 
               &z, &x, &y)
   return z
end

divexact{S, M}(x::PowerSeries{Residue{ZZ, M}, S}, y::Int) = divexact(x, ZZ(y))

/{T<: Ring, S}(x::PowerSeries{T, S}, y::Int) = divexact(x, y)

/{T<: Ring, S}(x::PowerSeries{T, S}, y::ZZ) = divexact(x, y)

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{T<: Ring, S}(a::PowerSeries{T, S})
   a == 0 && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   a1 = coeff(a, 0)
   if a.prec == nothing
      a.length != 1 && error("Unable to invert infinite precision power series")
      return PowerSeries(PowerSeries{T, S}, [divexact(T(1), a1)], nothing)
   end
   d = Array(T, a.prec)
   if a.prec != 0
      d[1] = divexact(T(1), a1)
   end
   a1 = -a1
   for n = 2:a.prec
      s = coeff(a, 1)*d[n - 1]
      for i = 2:min(n, a.length) - 1
         s += coeff(a, i)*d[n - i]
      end
      d[n] = divexact(s, a1)
   end
   ainv = PowerSeries(PowerSeries{T, S}, d, a.prec)
   ainv.length = normalise(ainv, a.prec)
   return ainv
end

function inv{S}(a::PowerSeries{ZZ, S})
   a == 0 && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   if a.prec == nothing
      a1 = coeff(a, 0)
      length(a) != 1 && error("Unable to invert infinite precision power series")
      return PowerSeries(PowerSeries{ZZ, S}, [divexact(ZZ(1), a1)], nothing)
   end
   ainv = PowerSeries(PowerSeries{ZZ, S}, Array(ZZ, 0), a.prec)
   ccall((:fmpz_poly_inv_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &ainv, &a, a.prec)
   return ainv
end

function inv{S, M}(a::PowerSeries{Residue{ZZ, M}, S})
   a == 0 && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   if a.prec == nothing
      a1 = coeff(a, 0)
      length(a) != 1 && error("Unable to invert infinite precision power series")
      return PowerSeries(PowerSeries{Residue{ZZ, M}, S}, [inv(a1)], nothing)
   end
   ainv = PowerSeries(PowerSeries{Residue{ZZ, M}, S}, Array(Residue{ZZ, M}, 0), a.prec)
   ccall((:fmpz_mod_poly_inv_series, :libflint), Void, 
                (Ptr{PowerSeries}, Ptr{PowerSeries}, Int), 
               &ainv, &a, a.prec)
   return ainv
end

###########################################################################################
#
#   Special functions
#
###########################################################################################

function exp{T<: Ring, S}(a::PowerSeries{T, S})
   if a == 0
      return PowerSeries(PowerSeries{T, S}, [one(T)], a.prec)
   elseif a.prec == nothing
      error("Unable to compute exponential of infinite precision power series")
   end
   d = Array(T, a.prec)
   d[0 + 1] = exp(coeff(a, 0))
   len = length(a)
   for k = 1 : a.prec - 1
      s = zero(T)
      for j = 1 : min(k + 1, len) - 1
         s += j * coeff(a, j) * d[k - j + 1]
      end
      !isunit(T(k)) && error("Unable to divide in exp")
      d[k + 1] = divexact(s, k)
   end
   b = PowerSeries(PowerSeries{T, S}, d, a.prec)
   return b
end

###########################################################################################
#
#   PowerSeriesRing constructor
#
###########################################################################################

function PowerSeriesRing{T <: Ring}(::Type{T}, s::String)
   S = symbol(s)
   T1 = PowerSeries{T, S}
   T2 = T
   
   # Conversions and promotions

   eval(:(Base.convert(::Type{$T1}, x::$T) = PowerSeries($T1, [x], nothing)))
   eval(:(Base.promote_rule(::Type{$T1}, ::Type{$T}) = $T1))

   P = T2.parameters
   while length(P) > 0
      T2 = P[1]
      if isa(T2, DataType) && T2 <: Ring
         eval(:(Base.convert(::Type{$T1}, x::$T2) = PowerSeries($T1, [convert($T, x)], nothing)))
         eval(:(Base.promote_rule(::Type{$T1}, ::Type{$T2}) = $T1))
         P = T2.parameters
      else
         break
      end
   end

   eval(:(Base.convert(::Type{$T1}, x::Integer) = PowerSeries($T1, [convert($T, x)], nothing)))
   eval(:(Base.promote_rule{R <: Integer}(::Type{$T1}, ::Type{R}) = $T1))

   # (Type, gen) 

   return (PowerSeries{T, S}, PowerSeries(PowerSeries{T, S}, [T(0), T(1)], nothing))
end
