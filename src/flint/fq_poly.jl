################################################################################
#
#  fq_poly.jl: Flint fq_poly (Polynomials over FqFiniteField)
#
################################################################################

export fq_poly, FqPolyRing

export elem_type, base_ring, parent, var, check_parent, length, coeff, zero,
       one, gen, isgen, degree, deepcopy, canonical_unit, show, show_minus_one,
       +, -, *, ^, ==, truncate, mullow, reverse, shift_left, shift_right, div,
       rem, divrem, gcd, evaluate, compose, derivative, inflate, deflate,
       factor, factor_distinct_deg, fit!, setcoeff!, mul!, add!, sub!, addeq!,
       Polynomial


################################################################################
#
#  Type and parent object methods
#
################################################################################

elem_type(::FqPolyRing) = fq_poly

base_ring(a::FqPolyRing) = a.base_ring

parent(a::fq_poly) = a.parent

var(a::FqPolyRing) = a.S

function check_parent(a::fq_poly, b::fq_poly) 
   a.parent != b.parent &&
         error("Operations on distinct polynomial rings not supported")
end

################################################################################
#
#   Basic manipulation
#
################################################################################
   
length(x::fq_poly) = ccall((:fq_poly_length, :libflint), Int,
                                (Ptr{fq_poly},), &x)

function coeff(x::fq_poly, n::Int)
   n < 0 && throw(DomainError())
   F = (x.parent).base_ring
   temp = F(1)
   ccall((:fq_poly_get_coeff, :libflint), Void, 
         (Ptr{fq}, Ptr{fq_poly}, Clong, Ptr{FqFiniteField}),
         &temp, &x, n, &F)
   return temp
end

zero(a::FqPolyRing) = a(zero(base_ring(a)))

one(a::FqPolyRing) = a(one(base_ring(a)))

gen(a::FqPolyRing) = a([zero(base_ring(a)), one(base_ring(a))])

isgen(x::fq_poly) = ccall((:fq_poly_is_gen, :libflint), Bool,
                              (Ptr{fq_poly}, Ptr{FqFiniteField}),
                              &x, &base_ring(x.parent))

degree(f::fq_poly) = f.length - 1

function deepcopy(a::fq_poly)
   z = fq_poly(a)
   z.parent = a.parent
   return z
end

################################################################################
#
#   Canonicalisation
#
################################################################################

canonical_unit(a::fq_poly) = canonical_unit(lead(a))
  
################################################################################
#
#  String I/O
#
################################################################################

function show(io::IO, x::fq_poly)
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fq_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                  (Ptr{fq_poly}, Ptr{Uint8}, Ptr{FqFiniteField}),
                  &x, bytestring(string(var(parent(x)))),
                  &((x.parent).base_ring))
      print(io, bytestring(cstr))
      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
end

function show(io::IO, R::FqPolyRing)
  print(io, "Univariate Polynomial Ring in ")
  print(io, string(var(R)))
  print(io, " over ")
  show(io, base_ring(R))
end

show_minus_one(::Type{fq_poly}) = show_minus_one(fq)

################################################################################
#
#  Uary operations
#
################################################################################

function -(x::fq_poly)
   z = parent(x)()
   ccall((:fq_poly_neg, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{FqFiniteField}),
         &z, &x, &base_ring(parent(x)))
   return z
end

################################################################################
#
#  Binary operations
#
################################################################################

function +(x::fq_poly, y::fq_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_poly_add, :libflint), Void, 
         (Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{fq_poly}, Ptr{FqFiniteField}),
         &z, &x, &y, &base_ring(parent(x)))
   return z
end

function -(x::fq_poly, y::fq_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_poly_sub, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{fq_poly}, Ptr{FqFiniteField}),
         &z, &x, &y, &base_ring(parent(x)))
   return z
end

function *(x::fq_poly, y::fq_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_poly_mul, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{fq_poly}, Ptr{FqFiniteField}),
         &z, &x, &y, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Ad hoc binary operators
#
################################################################################

function *(x::fq, y::fq_poly)
   x.parent == base_ring(parent(y)) &&
         error("Coefficient rings must be equal")
   z = parent(y)()
   ccall((:fq_poly_scalar_mul_fq, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{fq}, Ptr{FqFiniteField}),
         &z, &y, &x, &base_ring(parent(x)))
  return z
end

*(x::fq_poly, y::fq) = y*x

*(x::fmpz, y::fq) = base_ring(parent(y))(x) * y

*(x::fq, y::fmpz) = y*x

*(x::Integer, y::fq) = ZZ(x)*y

*(x::fq, y::Integer) = y*x

+(x::fq, y::fq_poly) = parent(y)(x) + y

+(x::fq_poly, y::fq) = y + x

+(x::fmpz, y::fq_poly) = base_ring(parent(y))(x) + y

+(x::fq_poly, y::fmpz) = y + x

+(x::fq_poly, y::Integer) = x + ZZ(y)

+(x::Integer, y::fq_poly) = y + x

################################################################################
#
#   Powering
#
################################################################################

function ^(x::fq_poly, y::Clong)
   y < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_poly_pow, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Clong, Ptr{FqFiniteField}), 
         &z, &x, y, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Comparisons
#
################################################################################

function ==(x::fq_poly, y::fq_poly)
   check_parent(x,y)
   r = ccall((:fq_poly_equal, :libflint), Cint,
             (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{FqFiniteField}),
             &x, &y, &base_ring(parent(x)))
   return Bool(r)
end

################################################################################
#
#   Ad hoc comparisons
#
################################################################################

function ==(x::fq_poly, y::fq) 
   base_ring(parent(x)) != parent(y) && return false
   if length(x) > 1
      return false
   elseif length(x) == 1 
      r = ccall((:fq_poly_equal_fq, :libflint), Cint, 
                (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{FqFiniteField}),
                &z, &temp, &base_ring(parent(x)))
      return Bool(r)
   else
      return y == 0
  end 
end

==(x::fq, y::fq_poly) = y == x

==(x::fq_poly, y::fmpz) = x == base_ring(parent(x))(y)

==(x::fmpz, y::fq_poly) = y == x

==(x::fq_poly, y::Integer) = x == ZZ(y)

==(x::Integer, y::fq_poly) = y == x

################################################################################
#
#   Truncation
#
################################################################################

function truncate(x::fq_poly, n::Clong)
   n < 0 && throw(DomainError())
   if length(x) <= n
      return x
   end
   z = parent(x)()
   ccall((:fq_poly_set_trunc, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Int, Ptr{FqFiniteField}),
         &z, &x, n, &base_ring(parent(x)))
   return z
end

function mullow(x::fq_poly, y::fq_poly, n::Clong)
   check_parent(x,y)
   n < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_poly_mullow, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Clong, Ptr{FqFiniteField}),
         &z, &x, &y, n, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Reversal
#
################################################################################

function reverse(x::fq_poly, len::Clong)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_poly_reverse, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Clong, Ptr{FqFiniteField}),
         &z, &x, len, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Shifting
#
################################################################################

function shift_left(x::fq_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_poly_shift_left, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Clong, Ptr{FqFiniteField}),
         &z, &x, len, &base_ring(parent(x)))
   return z
end

function shift_right(x::fq_poly, len::Clong)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_poly_shift_right, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Clong, Ptr{FqFiniteField}),
         &z, &x, len, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Euclidean division
#
################################################################################

function div(x::fq_poly, y::fq_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_poly_div_basecase, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{FqFiniteField}), &z, &x, &y, &base_ring(parent(x)))
  return z
end

function rem(x::fq_poly, y::fq_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_poly_rem, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{FqFiniteField}), &z, &x, &y, &base_ring(parent(x)))
  return z
end

function divrem(x::fq_poly, y::fq_poly)
   check_parent(x,y)
   z = parent(x)()
   r = parent(x)()
   ccall((:fq_poly_divrem, :libflint), Void, (Ptr{fq_poly},
         Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{FqFiniteField}), &z, &r, &x, &y, &base_ring(parent(x)))
   return z,r
end

################################################################################
#
#   GCD
#
################################################################################

function gcd(x::fq_poly, y::fq_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_poly_gcd, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{FqFiniteField}), &z, &x, &y, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Evaluation
#
################################################################################

function evaluate(x::fq_poly, y::fq)
   base_ring(parent(x)) != parent(y) && error("Incompatible coefficient rings")
   z = parent(y)()
   ccall((:fq_poly_evaluate_fq, :libflint), Void,
         (Ptr{fq}, Ptr{fq_poly}, Ptr{fq},
         Ptr{FqFiniteField}), &z, &x, &y, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Composition
#
################################################################################

function compose(x::fq_poly, y::fq_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_poly_compose, :libflint), Void, 
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{FqFiniteField}), &z, &x, &y, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Derivative
#
################################################################################

function derivative(x::fq_poly)
   z = parent(x)()
   ccall((:fq_poly_derivative, :libflint), Void, 
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{FqFiniteField}),
         &z, &x, &base_ring(parent(x)))
   return z
end

################################################################################
#
#  Inflation and deflation
#
################################################################################

function inflate(x::fq_poly, n::Int)
   z = parent(x)()
   ccall((:fq_poly_inflate, :libflint), Void, (Ptr{fq_poly},
         Ptr{fq_poly}, Culong, Ptr{FqFiniteField}),
         &z, &x, UInt(n), &base_ring(parent(x)))
   return z
end

function deflate(x::fq_poly, n::Int)
   z = parent(x)()
   ccall((:fq_poly_deflate, :libflint), Void,
         (Ptr{fq_poly}, Ptr{fq_poly}, Culong, Ptr{FqFiniteField}),
         &z, &x, UInt(n), &base_ring(parent(x)))
  return z
end

################################################################################
#
#  Factorization
#
################################################################################

function factor(x::fq_poly)
   R = parent(x)
   F = base_ring(R)
   a = F()
   fac = fq_poly_factor(F)
   ccall((:fq_poly_factor, :libflint), Void, (Ptr{fq_poly_factor},
         Ptr{fq}, Ptr{fq_poly}, Ptr{FqFiniteField}),
         &fac, &a, &x, &F)
   res = Array(Tuple{fq_poly,Clong},fac.num)
   for i in 1:fac.num
      f = R()
      ccall((:fq_poly_factor_get_poly, :libflint), Void,
            (Ptr{fq_poly}, Ptr{fq_poly_factor}, Clong,
            Ptr{FqFiniteField}), &f, &fac, i-1, &F)
      e = unsafe_load(fac.exp,i)
      res[i] = (f,e)
   end
   return res 
end  

function factor_distinct_deg(x::fq_poly)
   R = parent(x)
   F = base_ring(R)
   fac = fq_poly_factor(F)
   tmp = fq_poly_factor(F)
   ccall((:fq_poly_factor_fit_length, :libflint), Void,
         (Ptr{fq_poly_factor}, Clong, Ptr{FqFiniteField}),
         &tmp, degree(x), &F)
   ccall((:fq_poly_factor_distinct_deg, :libflint), Void, 
         (Ptr{fq_poly_factor}, Ptr{fq_poly}, Ptr{Clong},
         Ptr{FqFiniteField}), &fac, &x, &tmp.exp, &F)
   res = Array(Tuple{fq_poly, Clong}, fac.num)
   for i in 1:fac.num
      f = R()
      ccall((:fq_poly_factor_get_poly, :libflint), Void,
            (Ptr{fq_poly}, Ptr{fq_poly_factor}, Clong,
            Ptr{FqFiniteField}), &f, &fac, i-1, &F)
      d = unsafe_load(tmp.exp,i)
      res[i] = (f,d)
   end
   return res
end

################################################################################
#
#   Unsafe functions
#
################################################################################

function fit!(z::fq_poly, n::Clong)
   ccall((:fq_poly_fit_length, :libflint), Void, 
         (Ptr{fq_poly}, Clong, Ptr{FqFiniteField}),
         &z, n, &base_ring(parent(x)))
end

function setcoeff!(z::fq_poly, n::Int, x::fq)
   ccall((:fq_poly_set_coeff_fmpz, :libflint), Void, 
         (Ptr{fq_poly}, Int, Ptr{fq}, Ptr{FqFiniteField}),
         &z, n, &temp, &base_ring(parent(x)))
end

function mul!(z::fq_poly, x::fq_poly, y::fq_poly)
   ccall((:fq_poly_mul, :libflint), Void, 
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{FqFiniteField}), &z, &x, &y, &base_ring(parent(x)))
end

function add!(z::fq_poly, x::fq_poly, y::fq_poly)
   ccall((:fq_poly_add, :libflint), Void, 
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{FqFiniteField}), &z, &x, &y, &base_ring(parent(x)))
end

function sub!(z::fq_poly, x::fq_poly, y::fq_poly)
   ccall((:fq_poly_sub, :libflint), Void, 
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{FqFiniteField}), &z, &x, &y, &base_ring(parent(x)))
end


function addeq!(z::fq_poly, x::fq_poly)
   ccall((:fq_poly_add, :libflint), Void, 
         (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{fq_poly},
         Ptr{FqFiniteField}), &z, &z, &x, &base_ring(parent(x)))
end

################################################################################
#
#   Parent object call overloads
#
################################################################################

function Base.call(R::FqPolyRing)
   z = fq_poly()
   z.parent = R
   return z
end

function Base.call(R::FqPolyRing, x::fq)
  z = fq_poly(x)
  z.parent = R
  return z
end

function Base.call(R::FqPolyRing, x::fmpz)
   return R(base_ring(R)(x))
end

function Base.call(R::FqPolyRing, x::Integer)
   return R(ZZ(x))
end

function Base.call(R::FqPolyRing, x::Array{fq, 1})
   length(x) == 0 && error("Array must be non-empty")
   base_ring(R) != parent(x[1]) && error("Coefficient rings must coincide")
   z = fq_poly(x)
   z.parent = R
   return z
end

function Base.call(R::FqPolyRing, x::Array{fmpz, 1})
   length(x) == 0 && error("Array must be non-empty")
   z = fq_poly(x, base_ring(R))
   z.parent = R
   return z
end

function Base.call{T <: Integer}(R::FqPolyRing, x::Array{T, 1})
   length(x) == 0 && error("Array must be non-empty")
   return R(map(ZZ, x))
end

function Base.call(R::FqPolyRing, x::fmpz_poly)
   z = fq_poly(x, base_ring(R))
   z.parent = R
   return z
end

################################################################################
#
#   PolynomialRing constructor
#
################################################################################

function PolynomialRing(R::FqFiniteField, s::String)
   S = symbol(s)
   parent_obj = FqPolyRing(R, S)
   return parent_obj, parent_obj([R(0), R(1)])
end

