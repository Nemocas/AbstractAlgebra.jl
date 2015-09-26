###############################################################################
#
#   padic.jl : flint padic numbers
#
###############################################################################

export FlintPadicField, padic, prime, teichmuller, log, sqrt

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

function O(R::FlintPadicField, n::fmpz)
   if n == 1
      N = 0
   else
      p = prime(R)
      if n == p
         N = 1
      else
         N = flog(n, p) 
         p^(N) != n && error("Not a power of p in p-adic O()")
      end
   end
   d = padic(N)
   d.parent = R
   return d
end

function O(R::FlintPadicField, n::fmpq)
   m = den(n)
   if m == 1
      return O(R, num(n))
   end
   num(n) != 1 && error("Not a power of p in p-adic O()")
   p = prime(R)
   if m == p
      N = -1
   else
     N = -flog(m, p) 
     p^(-N) != m && error("Not a power of p in p-adic O()")
   end
   d = padic(N)
   d.parent = R
   return d
end

O(R::FlintPadicField, n::Integer) = O(R, fmpz(n))

elem_type(::FlintPadicField) = padic

base_ring(a::FlintPadicField) = Union{}

base_ring(a::padic) = Union{}

parent(a::padic) = a.parent

function check_parent(a::padic, b::padic)
   parent(a) != parent(b) && 
      error("Incompatible padic rings in padic operation")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function prime(R::FlintPadicField)
   z = fmpz()
   ccall((:padic_ctx_pow_ui, :libflint), Void, 
         (Ptr{fmpz}, Int, Ptr{FlintPadicField}), &z, 1, &R)
   return z 
end

precision(a::padic) = a.N

valuation(a::padic) = a.v

function zero(R::FlintPadicField)
   z = padic(R.prec_max)
   ccall((:padic_zero, :libflint), Void, (Ptr{padic},), &z)
   z.parent = R
   return z
end

function one(R::FlintPadicField)
   z = padic(R.prec_max)
   ccall((:padic_one, :libflint), Void, (Ptr{padic},), &z)
   z.parent = R
   return z
end

iszero(a::padic) = Bool(ccall((:padic_is_zero, :libflint), Cint, 
                              (Ptr{padic},), &a))

isone(a::padic) = Bool(ccall((:padic_is_one, :libflint), Cint,
                             (Ptr{padic},), &a))

###############################################################################
#
#   AbstractString{} I/O
#
###############################################################################

function show(io::IO, x::padic)
   ctx = parent(x)
   cstr = ccall((:padic_get_str, :libflint), Ptr{UInt8}, 
               (Ptr{Void}, Ptr{padic}, Ptr{FlintPadicField}),
                   C_NULL, &x, &ctx)

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{UInt8},), cstr)
   print(io, " + O(")
   print(io, prime(ctx))
   print(io, "^$(x.N))")
end

function show(io::IO, R::FlintPadicField)
   print(io, "Field of ", prime(R), "-adic numbers")
end

needs_parentheses(x::padic) = true

is_negative(x::padic) = false

show_minus_one(::FlintPadicField) = true

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::padic) = x

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(x::padic)
   if iszero(x)
      return x
   end
   ctx = parent(x)
   z = padic(x.N)
   ccall((:padic_neg, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), 
                     &z, &x, &ctx)
   z.parent = ctx
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(x::padic, y::padic)
   check_parent(x, y)
   ctx = parent(x)
   z = padic(min(x.N, y.N))
   z.parent = ctx
   ccall((:padic_add, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), 
               &z, &x, &y, &ctx)
   return z
end

function -(x::padic, y::padic)
   check_parent(x, y)
   ctx = parent(x)
   z = padic(min(x.N, y.N))
   z.parent = ctx
   ccall((:padic_sub, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), 
                  &z, &x, &y, &ctx)
   return z
end

function *(x::padic, y::padic)
   check_parent(x, y)
   ctx = parent(x)
   z = padic(min(x.N + y.v, y.N + x.v))
   z.parent = ctx
   ccall((:padic_mul, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), 
               &z, &x, &y, &ctx)
   return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::padic, b::padic)
   check_parent(a, b)
   ctx = parent(a)
   z = padic(min(a.N, b.N))
   ccall((:padic_sub, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), 
               &z, &a, &b, &ctx)
   return Bool(ccall((:padic_is_zero, :libflint), Cint, 
                (Ptr{padic},), &z))
end

function isequal(a::padic, b::padic)
   if parent(a) != parent(b)
      return false
   end
   return a.N == b.N && a == b
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::padic, n::Int)
   ctx = parent(a)
   z = padic(a.N + (n - 1)*a.v)
   z.parent = ctx
   ccall((:padic_pow_si, :libflint), Void, 
                (Ptr{padic}, Ptr{padic}, Int, Ptr{FlintPadicField}), 
               &z, &a, n, &ctx)
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::padic, b::padic)
   b == 0 && throw(DivideError())
   check_parent(a, b)
   ctx = parent(a)
   z = padic(min(a.N - b.v, b.N - 2*b.v + a.v))
   z.parent = ctx
   ccall((:padic_div, :libflint), Cint, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}),
               &z, &a, &b, &ctx)
   return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

divexact(a::padic, b::Int) = a*(fmpz(1)//fmpz(b))

divexact(a::padic, b::fmpz) = a*(1//b)

divexact(a::padic, b::fmpq) = a*inv(b)

divexact(a::Int, b::padic) = fmpz(a)*inv(b)

divexact(a::fmpz, b::padic) = inv((fmpz(1)//a)*b)

divexact(a::fmpq, b::padic) = inv(inv(a)*b)

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::padic)
   a == 0 && throw(DivideError())
   ctx = parent(a)
   z = padic(a.N - 2*a.v)
   z.parent = ctx
   ccall((:padic_inv, :libflint), Cint, 
         (Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), &z, &a, &ctx)
   return z
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(x::padic, y::padic)
   check_parent(x, y)
   if x == 0 && y == 0 
      z = zero(parent(x))
   else
      z = one(parent(x))
   end
   return z
end

###############################################################################
#
#   Square root
#
###############################################################################

function sqrt(a::padic)
   (a.v % 2) != 0 && error("Unable to take padic square root")
   ctx = parent(a)
   z = padic(a.N - div(a.v, 2))
   z.parent = ctx
   res = Bool(ccall((:padic_sqrt, :libflint), Cint, 
                    (Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), &z, &a, &ctx))      
   !res && error("Square root of p-adic does not exist")
   return z
end

###############################################################################
#
#   Special functions
#
###############################################################################

function exp(a::padic) 
   a != 0 && a.v <= 0 && throw(DomainError())
   ctx = parent(a)
   z = padic(a.N)
   z.parent = ctx
   res = Bool(ccall((:padic_exp, :libflint), Cint, 
                    (Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), &z, &a, &ctx))
   !res && error("Unable to compute exponential")
   return z
end

function log(a::padic) 
   (a.v > 0 || a.v < 0 || a == 0) && throw(DomainError())
   ctx = parent(a)
   z = padic(a.N)
   z.parent = ctx
   res = Bool(ccall((:padic_log, :libflint), Cint,
                    (Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), &z, &a, &ctx))
   !res && error("Unable to compute logarithm")
   return z
end

function teichmuller(a::padic)
   a.v < 0 && throw(DomainError())
   ctx = parent(a)
   z = padic(a.N)
   z.parent = ctx
   ccall((:padic_teichmuller, :libflint), Void,
         (Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), &z, &a, &ctx)
   return z
end
  
###############################################################################
#
#   Unsafe operators
#
###############################################################################

function mul!(z::padic, x::padic, y::padic)
   z.N = min(x.N + y.v, y.N + x.v)
   ctx = parent(x)
   ccall((:padic_mul, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), 
               &z, &x, &y, &ctx)
end

function addeq!(x::padic, y::padic)
   x.N = min(x.N, y.N)
   ctx = parent(x)
   ccall((:padic_add, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{FlintPadicField}), 
               &x, &x, &y, &ctx)
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

promote_rule{T<:Integer}(::Type{padic}, ::Type{T}) = padic

promote_rule(::Type{padic}, ::Type{fmpz}) = padic

promote_rule(::Type{padic}, ::Type{fmpq}) = padic

###############################################################################
#
#   Parent object overloads
#
###############################################################################

function Base.call(R::FlintPadicField)
   z = padic(R.prec_max)
   z.parent = R
   return z
end

function Base.call(R::FlintPadicField, n::fmpz)
   if n == 1
      N = 0
   else
      p = prime(R)
      N, = remove(n, p) 
   end
   z = padic(N + R.prec_max)
   ccall((:padic_set_fmpz, :libflint), Void, 
         (Ptr{padic}, Ptr{fmpz}, Ptr{FlintPadicField}), &z, &n, &R)
   z.parent = R
   return z
end

function Base.call(R::FlintPadicField, n::fmpq)
   m = den(n)
   if m == 1
      return R(num(n))
   end
   p = prime(R)
   if m == p
      N = -1
   else
     N = -flog(m, p) 
   end
   z = padic(N + R.prec_max)
   ccall((:padic_set_fmpq, :libflint), Void, 
         (Ptr{padic}, Ptr{fmpq}, Ptr{FlintPadicField}), &z, &n, &R)
   z.parent = R
   return z
end

Base.call(R::FlintPadicField, n::Integer) = R(fmpz(n))

###############################################################################
#
#   FlintPadicField constructor
#
###############################################################################

# inner constructor is also used directly

function FlintPadicField(p::Integer, prec::Int)
   return FlintPadicField(fmpz(p), prec)
end
