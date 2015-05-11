###############################################################################
#
#   fq_series.jl : Power series over flint finite fields
#
###############################################################################

export fq_series, FqSeriesRing

###############################################################################
#
#   Data types and memory management
#
###############################################################################

FqSeriesID = ObjectIdDict()

type FqSeriesRing <: Ring
   base_ring::Ring
   prec_max::Int
   S::Symbol

   function FqSeriesRing(R::Ring, prec::Int, s::Symbol)
      return try
         FqSeriesID[R, prec, s]
      catch
         FqSeriesID[R, prec, s] = new(R, prec, s)
      end
   end
end

type fq_series <: PowerSeriesElem
   coeffs::Ptr{Void}
   alloc::Int
   length::Int
   prec :: Int
   parent::FqSeriesRing

   function fq_series(ctx::FqFiniteField)
      z = new()
      ccall((:fq_poly_init, :libflint), Void, 
            (Ptr{fq_series}, Ptr{FqFiniteField}), &z, &ctx)
      finalizer(z, _fq_series_clear_fn)
      return z
   end
   
   function fq_series(ctx::FqFiniteField, a::Array{fq, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fq_poly_init2, :libflint), Void, 
            (Ptr{fq_series}, Int, Ptr{FqFiniteField}), &z, len, &ctx)
      for i = 1:len
         ccall((:fq_poly_set_coeff, :libflint), Void, 
               (Ptr{fq_series}, Int, Ptr{fq}, Ptr{FqFiniteField}),
                                               &z, i - 1, &a[i], &ctx)
      end
      z.prec = prec
      finalizer(z, _fq_series_clear_fn)
      return z
   end
   
   function fq_series(ctx::FqFiniteField, a::fq_series)
      z = new()
      ccall((:fq_poly_init, :libflint), Void, 
            (Ptr{fq_series}, Ptr{FqFiniteField}), &z, &ctx)
      ccall((:fq_poly_set, :libflint), Void, 
            (Ptr{fq_series}, Ptr{fq_series}, Ptr{FqFiniteField}), &z, &a, &ctx)
      finalizer(z, _fq_series_clear_fn)
      return z
   end
end

function _fq_series_clear_fn(a::fq_series)
   ctx = base_ring(a)
   ccall((:fq_poly_clear, :libflint), Void,
         (Ptr{fq_series}, Ptr{FqFiniteField}), &a, &ctx)
end

function O(a::fq_series)
   prec = length(a) - 1
   prec < 0 && throw(DomainError())
   z = fq_series(base_ring(a), Array(fq, 0), 0, prec)
   z.parent = parent(a)
   return z
end

elem_type(::FqSeriesRing) = fq_series

base_ring(R::FqSeriesRing) = R.base_ring

var(a::FqSeriesRing) = a.S

###############################################################################
#
#   Basic manipulation
#
###############################################################################    
   
max_precision(R::FqSeriesRing) = R.prec_max

function normalise(a::fq_series, len::Int)
   ctx = base_ring(a)
   if len > 0
      c = ctx()
      ccall((:fq_poly_get_coeff, :libflint), Void, 
         (Ptr{fq}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
                                           &c, &a, len - 1, &ctx)
   end
   while len > 0 && iszero(c)
      len -= 1
      if len > 0
         ccall((:fq_poly_get_coeff, :libflint), Void, 
            (Ptr{fq}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
                                           &c, &a, len - 1, &ctx)
      end
   end

   return len
end

function coeff(x::fq_series, n::Int)
   ctx = base_ring(x)
   if n < 0
      return ctx()
   end
   z = ctx()
   ccall((:fq_poly_get_coeff, :libflint), Void, 
         (Ptr{fq}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), &z, &x, n, &ctx)
   return z
end

zero(R::FqSeriesRing) = R(0)

one(R::FqSeriesRing) = R(1)

function gen(R::FqSeriesRing)
   ctx = base_ring(R)
   z = fq_series(ctx, [ctx(0), ctx(1)], 2, max_precision(R) + 1)
   z.parent = R
   return z
end

function deepcopy(a::fq_series)
   z = fq_series(base_ring(a), a)
   z.prec = a.prec
   z.parent = parent(a)
   return z
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, x::fq_series)
   if length(x) == 0
      print(io, "0")
   else
      ctx = base_ring(x)
      cstr = ccall((:fq_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
        (Ptr{fq_series}, Ptr{Uint8}, Ptr{FqFiniteField}), 
                     &x, bytestring(string(var(parent(x)))), &ctx)

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
   print(io, "+O(", string(var(parent(x))), "^", x.prec, ")")
end

function show(io::IO, a::FqSeriesRing)
   print(io, "Univariate power series ring in ", var(a), " over ")
   show(io, base_ring(a))
end

show_minus_one(::Type{fq_series}) = show_minus_one(fq)

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(x::fq_series)
   ctx = base_ring(x)
   z = parent(x)()
   ccall((:fq_poly_neg, :libflint), Void, 
                (Ptr{fq_series}, Ptr{fq_series}, Ptr{FqFiniteField}), 
               &z, &x, &ctx)
   z.prec = x.prec
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::fq_series, b::fq_series)
   check_parent(a, b)
   ctx = base_ring(a)
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = parent(a)()
   z.prec = prec
   ccall((:fq_poly_add_series, :libflint), Void, 
     (Ptr{fq_series}, Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
               &z, &a, &b, lenz, &ctx)
   return z
end

function -(a::fq_series, b::fq_series)
   check_parent(a, b)
   ctx = base_ring(a)
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = parent(a)()
   z.prec = prec
   ccall((:fq_poly_sub_series, :libflint), Void, 
     (Ptr{fq_series}, Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
               &z, &a, &b, lenz, &ctx)
   return z
end

function *(a::fq_series, b::fq_series)
   check_parent(a, b)
   ctx = base_ring(a)
   lena = length(a)
   lenb = length(b)
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   z = parent(a)()
   z.prec = prec
      
   if lena == 0 || lenb == 0
      return z
   end

   lenz = min(lena + lenb - 1, prec)

   ccall((:fq_poly_mullow, :libflint), Void, 
     (Ptr{fq_series}, Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}),
               &z, &a, &b, lenz, &ctx)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function setcoeff!(z::fq_series, n::Int, x::fq)
   ctx = base_ring(z)
   ccall((:fq_poly_set_coeff, :libflint), Void, 
                (Ptr{fq_series}, Int, Ptr{fq}, Ptr{FqFiniteField}), 
               &z, n, &x, &ctx)
end

function mul!(z::fq_series, a::fq_series, b::fq_series)
   ctx = base_ring(z)
   lena = length(a)
   lenb = length(b)
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   lenz = min(lena + lenb - 1, prec)
   if lenz < 0
      lenz = 0
   end

   z.prec = prec
   ccall((:fq_poly_mullow, :libflint), Void, 
     (Ptr{fq_series}, Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}),
               &z, &a, &b, lenz, &ctx)
end

function addeq!(a::fq_series, b::fq_series)
   ctx = base_ring(a)
   lena = length(a)
   lenb = length(b)
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   a.prec = prec
   ccall((:fq_poly_add_series, :libflint), Void, 
     (Ptr{fq_series}, Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}),
               &a, &a, &b, lenz, &ctx)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::fq, y::fq_series)
   ctx = base_ring(y)
   z = parent(y)()
   z.prec = y.prec
   ccall((:fq_poly_scalar_mul_fq, :libflint), Void, 
         (Ptr{fq_series}, Ptr{fq_series}, Ptr{fq}, Ptr{FqFiniteField}), 
               &z, &y, &x, &ctx)
   return z
end

###############################################################################
#
#   Shifting
#
###############################################################################

function shift_left(x::fq_series, len::Int)
   len < 0 && throw(DomainError())
   ctx = base_ring(x)
   xlen = length(x)
   z = parent(x)()
   z.prec = x.prec + len
   ccall((:fq_poly_shift_left, :libflint), Void, 
                (Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
               &z, &x, len, &ctx)
   return z
end

function shift_right(x::fq_series, len::Int)
   len < 0 && throw(DomainError())
   ctx = base_ring(x)
   xlen = length(x)
   z = parent(x)()
   if len >= xlen
      z.prec = max(0, x.prec - len)
   else
      z.prec = x.prec - len
      ccall((:fq_poly_shift_right, :libflint), Void, 
                (Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
               &z, &x, len, &ctx)
   end
   return z
end

###############################################################################
#
#   Truncation
#
###############################################################################

function truncate(x::fq_series, prec::Int)
   prec < 0 && throw(DomainError())
   if x.prec <= prec
      return x
   end
   ctx = base_ring(x)
   z = parent(x)()
   z.prec = prec
   ccall((:fq_poly_set_trunc, :libflint), Void, 
                (Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
               &z, &x, prec, &ctx)
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::fq_series, b::Int)
   b < 0 && throw(DomainError())
   ctx = base_ring(a)
   if isgen(a)
      z = parent(a)()
      setcoeff!(z, b, ctx(1))
      z.prec = a.prec + b - 1
   elseif length(a) == 0
      z = parent(a)()
      z.prec = a.prec + (b - 1)*valuation(a)
   elseif length(a) == 1
      return parent(a)([coeff(a, 0)^b], 1, a.prec)
   elseif b == 0
      return parent(a)([ctx(1)], 1, parent(a).prec_max)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
   end
   return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::fq_series, y::fq_series)
   check_parent(x, y)
   ctx = base_ring(x)
   prec = min(x.prec, y.prec)
   
   n = max(length(x), length(y))
   n = min(n, prec)
   
   return Bool(ccall((:fq_poly_equal_trunc, :libflint), Cint, 
                (Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
               &x, &y, n, &ctx))
end

function isequal(x::fq_series, y::fq_series)
   check_parent(x, y)
   ctx = base_ring(x)
   if x.prec != y.prec || length(x) != length(y)
      return false
   end
   return Bool(ccall((:fq_poly_equal, :libflint), Cint, 
                (Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
               &x, &y, length(x), &ctx))
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::fq_series, y::fq_series)
   check_parent(x, y)
   ctx = base_ring(x)
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
   z = parent(x)()
   z.prec = prec
   ccall((:fq_poly_div_series, :libflint), Void, 
     (Ptr{fq_series}, Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}),
               &z, &x, &y, prec, &ctx)
   return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::fq_series, y::fq)
   y == 0 && throw(DivideError())
   ctx = base_ring(x)
   z = parent(x)()
   z.prec = x.prec
   ccall((:fq_poly_scalar_div_fq, :libflint), Void, 
                (Ptr{fq_series}, Ptr{fq_series}, Ptr{fq}, Ptr{FqFiniteField}), 
               &z, &x, &y, &ctx)
   return z
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::fq_series)
   a == 0 && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   ctx = base_ring(a)
   ainv = parent(a)()
   ainv.prec = a.prec
   ccall((:fq_poly_inv_series, :libflint), Void, 
                (Ptr{fq_series}, Ptr{fq_series}, Int, Ptr{FqFiniteField}), 
               &ainv, &a, a.prec, &ctx)
   return ainv
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

Base.promote_rule{T <: Integer}(::Type{fq_series}, ::Type{T}) = fq_series

Base.promote_rule(::Type{fq_series}, ::Type{fmpz}) = fq_series

Base.promote_rule(::Type{fq_series}, ::Type{fq}) = fq_series

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call(a::FqSeriesRing)
   ctx = base_ring(a)
   z = fq_series(ctx)
   z.prec = a.prec_max
   z.parent = a
   return z
end

function Base.call(a::FqSeriesRing, b::Integer)
   ctx = base_ring(a)
   if b == 0
      z = fq_series(ctx)
      z.prec = a.prec_max
   else
      z = fq_series(ctx, [ctx(b)], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function Base.call(a::FqSeriesRing, b::fmpz)
   ctx = base_ring(a)
   if b == 0
      z = fq_series(ctx)
      z.prec = a.prec_max
   else
      z = fq_series(ctx, [ctx(b)], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function Base.call(a::FqSeriesRing, b::fq)
   ctx = base_ring(a)
   if b == 0
      z = fq_series(ctx)
      z.prec = a.prec_max
   else
      z = fq_series(ctx, [b], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function Base.call(a::FqSeriesRing, b::fq_series)
   parent(b) != a && error("Unable to coerce power series")
   return b
end

function Base.call(a::FqSeriesRing, b::Array{fq, 1}, len::Int, prec::Int)
   ctx = base_ring(a)
   z = fq_series(ctx, b, len, prec)
   z.parent = a
   return z
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

function PowerSeriesRing(R::FqFiniteField, prec::Int, s::String)
   S = symbol(s)

   parent_obj = FqSeriesRing(R, prec, S)

   return parent_obj, parent_obj([R(0), R(1)], 2, prec)
end

