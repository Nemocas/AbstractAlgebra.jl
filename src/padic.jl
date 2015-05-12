###############################################################################
#
#   padic.jl : flint padic numbers
#
###############################################################################

import Base: log, sqrt

export PadicField, Padic, prime, teichmuller, log, sqrt

###############################################################################
#
#   Data types and memory management
#
###############################################################################

PadicBase = ObjectIdDict()

type PadicField <: Field
   p::Int 
   pinv::Float64
   pow::Ptr{Void}
   minpre::Int
   maxpre::Int
   mode::Int
   prec_max::Int

   function PadicField(prec::Int, p::fmpz)
      !isprime(p) && error("Prime base required in PadicField")
      d = new()
      ccall((:padic_ctx_init, :libflint), Void, 
           (Ptr{PadicField}, Ptr{fmpz}, Int, Int, Cint), 
                                     &d, &p, 0, 0, 1)
      finalizer(d, _padic_ctx_clear_fn)
      d.prec_max = prec
      return d
   end
end

function _padic_ctx_clear_fn(a::PadicField)
   ccall((:padic_ctx_clear, :libflint), Void, (Ptr{PadicField},), &a)
end

type padic <: FieldElem
   u :: Int
   v :: Int
   N :: Int
   parent::PadicField

   function padic(prec::Int)
      d = new()
      ccall((:padic_init2, :libflint), Void, (Ptr{padic}, Int), &d, prec)
      finalizer(d, _padic_clear_fn)
      return d
   end
end

function _padic_clear_fn(a::padic)
   ccall((:padic_clear, :libflint), Void, (Ptr{padic},), &a)
end

function O(R::PadicField, n::fmpz)
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

function O(R::PadicField, n::fmpq)
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

O(R::PadicField, n::Integer) = O(R, ZZ(n))

elem_type(::PadicField) = padic

base_ring(a::PadicField) = None

base_ring(a::padic) = None

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

function prime(R::PadicField)
   z = ZZ()
   ccall((:padic_ctx_pow_ui, :libflint), Void, 
         (Ptr{fmpz}, Int, Ptr{PadicField}), &z, 1, &R)
   return z 
end

precision(a::padic) = a.N

valuation(a::padic) = a.v

function zero(R::PadicField)
   z = padic(R.prec_max)
   ccall((:padic_zero, :libflint), Void, (Ptr{padic},), &z)
   z.parent = R
   return z
end

function one(R::PadicField)
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
#   String I/O
#
###############################################################################

function show(io::IO, x::padic)
   ctx = parent(x)
   cstr = ccall((:padic_get_str, :libflint), Ptr{Uint8}, 
               (Ptr{Void}, Ptr{padic}, Ptr{PadicField}),
                   C_NULL, &x, &ctx)

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   print(io, " + O(")
   print(io, prime(ctx))
   print(io, "^$(x.N))")
end

function show(io::IO, R::PadicField)
   print(io, "Field of ", prime(R), "-adic numbers")
end

needs_parentheses(x::padic) = true

is_negative(x::padic) = false

show_minus_one(::PadicField) = true

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
         (Ptr{padic}, Ptr{padic}, Ptr{PadicField}), 
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
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{PadicField}), 
               &z, &x, &y, &ctx)
   return z
end

function -(x::padic, y::padic)
   check_parent(x, y)
   ctx = parent(x)
   z = padic(min(x.N, y.N))
   z.parent = ctx
   ccall((:padic_sub, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{PadicField}), 
                  &z, &x, &y, &ctx)
   return z
end

function *(x::padic, y::padic)
   check_parent(x, y)
   ctx = parent(x)
   z = padic(min(x.N + y.v, y.N + x.v))
   z.parent = ctx
   ccall((:padic_mul, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{PadicField}), 
               &z, &x, &y, &ctx)
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
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{PadicField}), 
               &z, &x, &y, &ctx)
end

function addeq!(x::padic, y::padic)
   x.N = min(x.N, y.N)
   ctx = parent(x)
   ccall((:padic_add, :libflint), Void, 
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{PadicField}), 
               &x, &x, &y, &ctx)
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
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{PadicField}), 
               &z, &a, &b, &ctx)
   return Bool(ccall((:padic_is_zero, :libflint), Cint, 
                (Ptr{padic},), &z))
end

function isequal(a::padic, b::padic)
   check_parent(a, b)
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
                (Ptr{padic}, Ptr{padic}, Int, Ptr{PadicField}), 
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
         (Ptr{padic}, Ptr{padic}, Ptr{padic}, Ptr{PadicField}),
               &z, &a, &b, &ctx)
   return z
end

divexact(a::padic, b::Int) = a*(ZZ(1)//ZZ(b))

divexact(a::padic, b::fmpz) = a*(1//b)

divexact(a::padic, b::fmpq) = a*inv(b)

divexact(a::Int, b::padic) = ZZ(a)*inv(b)

divexact(a::fmpz, b::padic) = inv((ZZ(1)//a)*b)

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
         (Ptr{padic}, Ptr{padic}, Ptr{PadicField}), &z, &a, &ctx)
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
                    (Ptr{padic}, Ptr{padic}, Ptr{PadicField}), &z, &a, &ctx))      
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
                    (Ptr{padic}, Ptr{padic}, Ptr{PadicField}), &z, &a, &ctx))
   !res && error("Unable to compute exponential")
   return z
end

function log(a::padic) 
   (a.v > 0 || a.v < 0 || a == 0) && throw(DomainError())
   ctx = parent(a)
   z = padic(a.N)
   z.parent = ctx
   res = Bool(ccall((:padic_log, :libflint), Cint,
                    (Ptr{padic}, Ptr{padic}, Ptr{PadicField}), &z, &a, &ctx))
   !res && error("Unable to compute logarithm")
   return z
end

function teichmuller(a::padic)
   a.v < 0 && throw(DomainError())
   ctx = parent(a)
   z = padic(a.N)
   z.parent = ctx
   ccall((:padic_teichmuller, :libflint), Void,
         (Ptr{padic}, Ptr{padic}, Ptr{PadicField}), &z, &a, &ctx)
   return z
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

function Base.call(R::PadicField)
   z = padic(R.prec_max)
   z.parent = R
   return z
end

function Base.call(R::PadicField, n::fmpz)
   if n == 1
      N = 0
   else
      p = prime(R)
      N, = remove(n, p) 
   end
   z = padic(N + R.prec_max)
   ccall((:padic_set_fmpz, :libflint), Void, 
         (Ptr{padic}, Ptr{fmpz}, Ptr{PadicField}), &z, &n, &R)
   z.parent = R
   return z
end

function Base.call(R::PadicField, n::fmpq)
   m = den(n)
   if m == 1
      return R(num(n))
   end
   p = prime(R)
   if m == p
      N = -1
   else
     N = -flog(m, p) 
     p^(-N) != m && error("Not a power of p in denominator of p-adic coercion")
   end
   z = padic(N + R.prec_max)
   ccall((:padic_set_fmpq, :libflint), Void, 
         (Ptr{padic}, Ptr{fmpq}, Ptr{PadicField}), &z, &n, &R)
   z.parent = R
   return z
end

Base.call(R::PadicField, n::Integer) = R(ZZ(n))

###############################################################################
#
#   PadicField constructor
#
###############################################################################

# inner constructor is also used directly

function PadicField(prec::Int, p::Integer)
   return PadicField(prec, ZZ(p))
end