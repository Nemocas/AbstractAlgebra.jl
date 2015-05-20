################################################################################
#
#  nmod_poly.jl : Flint nmod_poly
#
################################################################################

export NmodPolyRing, nmod_poly

export parent, base_ring, elem_type, length, zero, one, gen, isgen, iszero,
       var, deepcopy, show, -, +, *, ^, ==, truncate, mullow, reverse,
       shift_left, shift_right, divexact, divrem, rem, gcd, xgcd, resultant,
       evaluate, derivative, compose, interpolate, inflate, deflate, lift,
       isirreducible, issquarefree, factor, factor_squarefree,
       factor_distinct_deg, factor_shape, setcoeff!, canonical_unit,
       add!, sub!, mul!, call, PolynomialRing, check_parent, gcdx, mod,
       invmod, gcdinv, mulmod, powmod

################################################################################
#
#  Types and memory management
#
################################################################################

NmodPolyRingID = ObjectIdDict()

type NmodPolyRing <: Ring
  base_ring::ResidueRing{fmpz}
  S::Symbol
  _n::UInt

  function NmodPolyRing(R::ResidueRing{fmpz}, s::Symbol)
    ZZ(typemax(UInt)) < R.modulus &&
      error("Modulus of residue ring must be less then ", ZZ(typemax(UInt)))

    return try
       NmodPolyRingID[R, s]
    catch
       NmodPolyRingID[R, s] = new(R, s, UInt(R.modulus))
    end
  end
end

type nmod_poly <: PolyElem
   _coeffs::Ptr{Void}
   _alloc::Int
   _length::Int
   _mod_n::UInt64
   _mod_ninv::UInt64
   _mod_norm::UInt64
   parent::NmodPolyRing

   function nmod_poly(n::UInt64)
      z = new()
      ccall((:nmod_poly_init, :libflint), Void, (Ptr{nmod_poly}, UInt64), &z, n)
      finalizer(z, _nmod_poly_clear_fn)
      return z
   end

  function nmod_poly{T <: Integer}(n::UInt64, arr::Array{T, 1})
    length(arr) == 0 && error("Array must have length > 0")
    arr = map(ZZ, arr)
    return nmod_poly(n, arr)
  end

  function nmod_poly(n::UInt64, arr::Array{fmpz, 1})
    length(arr) == 0  && error("Array must have length > 0")
    z = new()
    t = ZZ()
    tt = UInt(0)
    ccall((:nmod_poly_init2, :libflint), Void,
            (Ptr{nmod_poly}, UInt64, Int), &z, n, length(arr))
    for i in 1:length(arr)
      ccall((:fmpz_mod_ui, :libflint), UInt,
              (Ptr{fmpz}, Ptr{fmpz}, UInt), &t, &arr[i], n)
      tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
      ccall((:nmod_poly_set_coeff_ui, :libflint), Void,
              (Ptr{nmod_poly}, Int, UInt64), &z, i-1, tt)
    end
    finalizer(z, _nmod_poly_clear_fn)
    return z
  end

  function nmod_poly(n::UInt64, arr::Array{Residue{fmpz}, 1})
    length(arr) == 0  && error("Array must have length > 0")
    ZZ(n) != modulus(arr[1]) && error("Moduli must coincide")
    z = new()
    t = ZZ()
    tt = UInt(0)
    ccall((:nmod_poly_init2, :libflint), Void,
            (Ptr{nmod_poly}, UInt64, Int), &z, n, length(arr))
    for i in 1:length(arr)
      ccall((:fmpz_mod_ui, :libflint), UInt,
                (Ptr{fmpz}, Ptr{fmpz}, UInt), &t, &arr[i].data, n)
      tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
      ccall((:nmod_poly_set_coeff_ui, :libflint), Void,
              (Ptr{nmod_poly}, Int, UInt64), &z, i-1, tt)
    end
    finalizer(z, _nmod_poly_clear_fn)
    return z
  end

  function nmod_poly(n::UInt64, f::fmpz_poly)
    z = new()
    ccall((:nmod_poly_init2, :libflint), Void,
            (Ptr{nmod_poly}, UInt64, Int), &z, n, degree(f))
    ccall((:fmpz_poly_get_nmod_poly, :libflint), Void,
            (Ptr{nmod_poly}, Ptr{fmpz_poly}), &z, &f)
    finalizer(z, _nmod_poly_clear_fn)
    return z
  end

  function nmod_poly(f::nmod_poly)
    z = new()
    ccall((:nmod_poly_init, :libflint), Void, (Ptr{nmod_poly}, ), &z)
    ccall((:nmod_poly_set, :libflint), Void,
            (Ptr{nmod_poly}, Ptr{nmod_poly}), &z, &f)
    finalizer(z, _nmod_poly_clear_fn)
    return z
  end
end

function _nmod_poly_clear_fn(x::nmod_poly)
  ccall((:nmod_poly_clear, :libflint), Void, (Ptr{nmod_poly}, ), &x)
end

type nmod_poly_factor
  poly::Ptr{nmod_poly}
  exp::Ptr{UInt64} # this is mp_limb_signed_t in gmp and slong in flint
  _num::Int
  _alloc::Int
  _n::UInt64
    
  function nmod_poly_factor(n::UInt64)
    z = new()
    ccall((:nmod_poly_factor_init, :libflint), Void,
            (Ptr{nmod_poly_factor}, ), &z)
    z._n = n
    finalizer(z, _nmod_poly_factor_clear_fn)
    return z
  end
end

function _nmod_poly_factor_clear_fn(a::nmod_poly_factor)
  ccall((:nmod_poly_factor_clear, :libflint), Void,
          (Ptr{nmod_poly_factor}, ), &a)
end

parent(a::nmod_poly) = a.parent

base_ring(R::NmodPolyRing) = R.base_ring

base_ring(a::nmod_poly) = base_ring(parent(a))

elem_type(::nmod_poly) = nmod_poly

elem_type(::NmodPolyRing) = nmod_poly

################################################################################
#
#  Basic manipulation
#
################################################################################

length(x::nmod_poly) = ccall((:nmod_poly_length, :libflint), Int,
                               (Ptr{nmod_poly}, ), &x)

degree(x::nmod_poly) = ccall((:nmod_poly_degree, :libflint), Int,
                               (Ptr{nmod_poly}, ), &x)

function coeff(x::nmod_poly, n::Int)
  (n < 0 || n > degree(x)) && throw(DomainError())
  return base_ring(x)(ccall((:nmod_poly_get_coeff_ui, :libflint), UInt64,
          (Ptr{nmod_poly}, Int), &x, n))
end

zero(R::NmodPolyRing) = R(0)

one(R::NmodPolyRing) = R([1])

gen(R::NmodPolyRing) = R([zero(base_ring(R)), one(base_ring(R))])

isgen(a::nmod_poly) = (degree(a) == 1 &&
                              iszero(coeff(a,0)) && isone(coeff(a,1)))

iszero(a::nmod_poly) = Bool(ccall((:nmod_poly_is_zero, :libflint), Int32,
                              (Ptr{nmod_poly}, ), &a))

var(R::NmodPolyRing) = R.S

function deepcopy(a::nmod_poly)
  z = nmod_poly(a)
  if isdefined(a, :parent)
    z.parent = a.parent
  end
  return z
end

################################################################################
#
#  String I/O
#
################################################################################

function show(io::IO, x::nmod_poly)
  if length(x) == 0
    print(io, "0")
  else
    cstr = ccall((:nmod_poly_get_str_pretty, :libflint), Ptr{Uint8},
            (Ptr{nmod_poly}, Ptr{Uint8}), &x, bytestring(string(var(parent(x)))))
    print(io, bytestring(cstr))
    ccall((:flint_free, :libflint), Void, (Ptr{Uint8}, ), cstr)
  end
end

function show(io::IO, R::NmodPolyRing)
  print(io, "Univariate Polynomial Ring in ")
  print(io, string(var(R)))
  print(io, " over ")
  print(io, base_ring(R))
end

################################################################################
#
#  Unary operations
#
################################################################################

function -(x::nmod_poly)
  z = parent(x)()
  ccall((:nmod_poly_neg, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}), &z, &x)
  return z
end

################################################################################
#
#   Binary operations
#
################################################################################

function +(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  z = parent(x)()
  ccall((:nmod_poly_add, :libflint), Void, 
          (Ptr{nmod_poly}, Ptr{nmod_poly},  Ptr{nmod_poly}), &z, &x, &y)
  return z
end

function -(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  z = parent(x)()
  ccall((:nmod_poly_sub, :libflint), Void, 
          (Ptr{nmod_poly}, Ptr{nmod_poly},  Ptr{nmod_poly}), &z, &x, &y)
  return z
end

function *(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  z = parent(x)()
  ccall((:nmod_poly_mul, :libflint), Void, 
          (Ptr{nmod_poly}, Ptr{nmod_poly},  Ptr{nmod_poly}), &z, &x, &y)
  return z
end

###############################################################################
#
#  Ad hoc binary operations
#
###############################################################################

function *(x::nmod_poly, y::UInt64)
  z = parent(x)()
  ccall((:nmod_poly_scalar_mul_nmod, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, UInt64), &z, &x, y)
  return z
end

*(x::UInt64, y::nmod_poly) = y*x

function *(x::nmod_poly, y::fmpz)
  z = parent(x)()
  t = ZZ()
  tt = UInt(0)
  ccall((:fmpz_mod_ui, :libflint), UInt,
                (Ptr{fmpz}, Ptr{fmpz}, UInt), &t, &y, parent(x)._n)
  tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
  return x*tt
end

*(x::fmpz, y::nmod_poly) = y*x

*(x::nmod_poly, y::Integer) = x*ZZ(y)

*(x::Integer, y::nmod_poly) = y*x

function *(x::nmod_poly, y::Residue{fmpz})
  (base_ring(x) != parent(y)) && error("Must have same parent")
  return x*y.data
end

*(x::Residue{fmpz}, y::nmod_poly) = y*x

function +(x::nmod_poly, y::UInt64)
  z = parent(x)()
  ccall((:nmod_poly_add_ui, :libflint), Void,
    (Ptr{nmod_poly}, Ptr{nmod_poly}, UInt64), &z, &x, y)
  return z
end

+(x::UInt64, y::nmod_poly) = +(y,x)

function +(x::nmod_poly, y::fmpz)
  z = parent(x)()
  t = ZZ()
  tt = UInt(0)
  ccall((:fmpz_mod_ui, :libflint), UInt,
                (Ptr{fmpz}, Ptr{fmpz}, UInt), &t, &y, parent(x)._n)
  tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
  return +(x,tt)
end

+(x::fmpz, y::nmod_poly) = +(y,x)

+(x::nmod_poly, y::Integer) = +(x,ZZ(y))

+(x::Integer, y::nmod_poly) = +(y,x) 

function +(x::nmod_poly, y::Residue{fmpz})
  (base_ring(x) != parent(y)) && error("Elements must have same parent")
  return +(x,y.data)
end

+(x::Residue{fmpz}, y::nmod_poly) = y + x

function -(x::nmod_poly, y::UInt64)
  z = parent(x)()
  ccall((:nmod_poly_sub_ui, :libflint), Void,
    (Ptr{nmod_poly}, Ptr{nmod_poly}, UInt64), &z, &x, y)
  return z
end

-(x::UInt64, y::nmod_poly) = -(-(y,x))

function -(x::nmod_poly, y::fmpz)
  z = parent(x)()
  t = ZZ()
  tt = UInt(0)
  ccall((:fmpz_mod_ui, :libflint), UInt,
                (Ptr{fmpz}, Ptr{fmpz}, UInt), &t, &y, parent(x)._n)
  tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
  return -(x,tt)
end

-(x::fmpz, y::nmod_poly) = -(-(y,x))

-(x::nmod_poly, y::Integer) = -(x,ZZ(y))

-(x::Integer, y::nmod_poly) = -(-(y,x))

function -(x::nmod_poly, y::Residue{fmpz})
  (base_ring(x) != parent(y)) && error("Elements must have same parent")
  return -(x,y.data)
end

-(x::Residue{fmpz}, y::nmod_poly) = -(-(y,x))

################################################################################
#
#  Powering
#
################################################################################

function ^(x::nmod_poly, y::UInt64)
  z = parent(x)()
  ccall((:nmod_poly_pow, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, UInt64), &z, &x, y)
  return z
end

function ^(x::nmod_poly, y::Int)
  y < 0 && throw(DomainError())
  return x^UInt(y)
end

function ^(x::nmod_poly, y::fmpz)
  y < 0 && throw(DomainError())
  y > ZZ(typemax(UInt)) &&
          error("Exponent must be smaller than ", typemax(UInt))
  return x^UInt(y)
end
 
################################################################################
#
#  Comparison
#
################################################################################

function ==(x::nmod_poly, y::nmod_poly)
  return parent(x) == parent(y) &&
  Bool(ccall((:nmod_poly_equal, :libflint), Int32,
          (Ptr{nmod_poly}, Ptr{nmod_poly}), &x, &y))
end

==(x::nmod_poly, y::Int) = x == parent(x)(y)

==(x::Int, y::nmod_poly) = y == x

################################################################################
#
#  Ad hoc comparisons
#
################################################################################

function ==(x::nmod_poly, y::Residue{fmpz})
  if base_ring(x) != parent(y)
    return false
  end
  if length(x) > 1
    return false
  elseif length(x) == 1 
    u = ccall((:nmod_poly_get_coeff_ui, :libflint), UInt64, 
            (Ptr{nmod_poly}, Int), &x, 0)
    return parent(y)(u) == y
  else
    return y == 0
  end 
end

==(x::Residue{fmpz}, y::nmod_poly) = y == x

################################################################################
#
#  Truncation
#
################################################################################

function truncate(a::nmod_poly, n::Int)
  n < 0 && throw(DomainError())

  z = deepcopy(a)
   
  if length(z) <= n
    return z
  end

  ccall((:nmod_poly_truncate, :libflint), Void,
          (Ptr{nmod_poly}, Int), &z, n)
  return z
end

function mullow(x::nmod_poly, y::nmod_poly, n::Int)
  check_parent(x, y)
  n < 0 && throw(DomainError())

  z = parent(x)()
  ccall((:nmod_poly_mullow, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}, Int), &z, &x, &y, n)
  return z
end

###############################################################################
#
#   Reversal
#
###############################################################################

function reverse(x::nmod_poly, len::Int)
  len < 0 && throw(DomainError())
  z = parent(x)()
  ccall((:nmod_poly_reverse, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Int), &z, &x, len)
  return z
end

###############################################################################
#
#   Shifting
#
###############################################################################

function shift_left(x::nmod_poly, len::Int)
  len < 0 && throw(DomainError())
  z = parent(x)()
  ccall((:nmod_poly_shift_left, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Int), &z, &x, len)
  return z
end

function shift_right(x::nmod_poly, len::Int)
  len < 0 && throw(DomainError())
  z = parent(x)()
  ccall((:nmod_poly_shift_right, :libflint), Void,
            (Ptr{nmod_poly}, Ptr{nmod_poly}, Int), &z, &x, len)
  return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::nmod_poly, y::nmod_poly)
  check_parent(x, y)
  iszero(y) && throw(DivideError())
  z = parent(x)()
  ccall((:nmod_poly_div, :libflint), Void, 
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}), &z, &x, &y)
  return z
end

################################################################################
#
#  Ad hoc exact division
#
################################################################################

function divexact(x::nmod_poly, y::Residue{fmpz})
  base_ring(x) != parent(y) && error("Elements must have same parent")
  iszero(y) && throw(DivideError())
  return divexact(x,parent(x)(y))
end

################################################################################
#
#  Division with remainder
#
################################################################################

function divrem(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  iszero(y) && throw(DivideError()) 
  q = parent(x)()
  r = parent(x)()
  ccall((:nmod_poly_divrem, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}),
          &q, &r, &x, &y)
  return q,r
end

################################################################################
#
#  Remainder
#
################################################################################

function rem(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  iszero(y) && throw(DivideError()) 
  z = parent(x)()
  ccall((:nmod_poly_rem, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}), &z, &x, &y)
  return z
end

mod(x::nmod_poly, y::nmod_poly) = rem(x,y)

################################################################################
#
#  GCD 
#
################################################################################

function gcd(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  z = parent(x)()
  ccall((:nmod_poly_gcd, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}), &z, &x, &y)
  return z
end 

function xgcd(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  g = parent(x)()
  s = parent(x)()
  t = parent(x)()
  ccall((:nmod_poly_xgcd, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly},
           Ptr{nmod_poly}), &g, &s, &t, &x, &y)
  return g,s,t
end

gcdx(x::nmod_poly, y::nmod_poly) = xgcd(x,y)

function gcdinv(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  length(y) <= 1 && error("Length of second argument must be >= 2")
  g = parent(x)()
  s = parent(x)()
  ccall((:nmod_poly_gcdinv, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}),
          &g, &s, &x, &y)
  return g,s
end 

################################################################################
#
#  Modular arithmetic
#
################################################################################

function invmod(x::nmod_poly, y::nmod_poly)
  length(y) == 0 && error("Second argument must not be 0")

  check_parent(x,y)

  if length(y) == 1 
    return parent(x)(inv(eval(x,coeff(y,0))))
  end

  z = parent(x)()

  r = ccall((:nmod_poly_invmod, :libflint), Cint,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}), &z, &x, &y)

  r == 0 ? error("GCD not 1") : return z
end

function mulmod(x::nmod_poly, y::nmod_poly, z::nmod_poly)
  check_parent(x,y)
  check_parent(y,z)

  w = parent(x)()

  ccall((:nmod_poly_mulmod, :libflint), Void,
        (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}),
        &w, &x, &y, &z)
  
  return w
end

function powmod(x::nmod_poly, e::UInt, y::nmod_poly)
  z = parent(x)()

  ccall((:nmod_poly_powmod_ui_binexp, :libflint), Void,
  (Ptr{nmod_poly}, Ptr{nmod_poly}, UInt, Ptr{nmod_poly}), &z, &x, e, &y)

  return z
end

function powmod(x::nmod_poly, e::Int, y::nmod_poly)
  e < 0 && error("Exponent must be positive")

  return powmod(x, UInt(e), y)
end

################################################################################
#
#  Resultant
#
################################################################################

function resultant(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  z = parent(x)()
  r = ccall((:nmod_poly_resultant, :libflint), UInt,
          (Ptr{nmod_poly}, Ptr{nmod_poly}), &x, &y)
  return base_ring(x)(r)
end

################################################################################
#
#  Evaluation
#
################################################################################

function evaluate(x::nmod_poly, y::Residue{fmpz})
  base_ring(x) != parent(y) && error("Elements must have same parent")
  t = ZZ()
  ccall((:fmpz_mod_ui, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, UInt64), &t, &y.data, parent(x)._n)
  u = ccall((:fmpz_get_ui, :libflint), UInt64, (Ptr{fmpz}, ), &t)
  z = ccall((:nmod_poly_evaluate_nmod, :libflint), UInt,
              (Ptr{nmod_poly}, UInt), &x, u)
  return parent(y)(z)
end

################################################################################
#
#  Derivative
#
################################################################################

function derivative(x::nmod_poly)
  z = parent(x)()
  ccall((:nmod_poly_derivative, :libflint), Void,
        (Ptr{nmod_poly}, Ptr{nmod_poly}), &z, &x)
  return z
end

################################################################################
#
#  Integral
#
################################################################################

function integral(x::nmod_poly)
  z = parent(x)()
  ccall((:nmod_poly_integral, :libflint), Void,
        (Ptr{nmod_poly}, Ptr{nmod_poly}), &z, &x)
  return z
end

################################################################################
#
#  Composition
#
################################################################################

function compose(x::nmod_poly, y::nmod_poly)
  check_parent(x,y)
  z = parent(x)()
  ccall((:nmod_poly_compose, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, Ptr{nmod_poly}), &z, &x, &y)
  return z
end

################################################################################
#
#  Interpolation
#
################################################################################

function interpolate(R::NmodPolyRing, x::Array{Residue{fmpz}, 1},
                                      y::Array{Residue{fmpz}, 1})
  z = R()

  ax = Array(UInt64, length(x))
  ay = Array(UInt64, length(y))

  t = ZZ()

  for i in 1:length(x)
    ccall((:fmpz_mod_ui, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, UInt64), &t, &x[i].data, R._n)
    u = ccall((:fmpz_get_ui, :libflint), UInt64, (Ptr{fmpz}, ), &t)
    ax[i] = u

    ccall((:fmpz_mod_ui, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, UInt64), &t, &y[i].data, R._n)
    u = ccall((:fmpz_get_ui, :libflint), UInt64, (Ptr{fmpz}, ), &t)
    
    ay[i] = u
  end
  ccall((:nmod_poly_interpolate_nmod_vec, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{UInt64}, Ptr{UInt64}, Int),
          &z, ax, ay, length(x))
  return z
end

################################################################################
#
#  Inflation and Deflation
#
################################################################################

function inflate(x::nmod_poly, n::Int)
  n < 0 && throw(DomainError())
  z = parent(x)()
  ccall((:nmod_poly_inflate, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, UInt64), &z, &x, UInt(n))
  return z
end

function deflate(x::nmod_poly, n::Int)
  n < 0 && throw(DomainError())
  z = parent(x)()
  ccall((:nmod_poly_deflate, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}, UInt64), &z, &x, UInt(n))
  return z
end
 
################################################################################
#
#  Lifting
#
################################################################################

function lift(x::FmpzPolyRing, y::nmod_poly)
  base_ring(x) != ZZ && error("Can only lift to integer polynomial ring")
  z = _lift(y)
  z.parent = x
  return z
end

################################################################################
#
#  Irreducibility
#
################################################################################

function isirreducible(x::nmod_poly)
  r = ccall((:nmod_poly_is_irreducible, :libflint), Int32,
          (Ptr{nmod_poly}, ), &x)
  return Bool(r)
end

################################################################################
#
#  Squarefree(ness)
#
################################################################################

function issquarefree(x::nmod_poly)
  r = ccall((:nmod_poly_is_squarefree, :libflint), Int32, (Ptr{nmod_poly}, ), &x)
  return Bool(r)
end

################################################################################
#
#  Factorization
#
################################################################################

function factor(x::nmod_poly)
  fac = nmod_poly_factor(x._mod_n)
  ccall((:nmod_poly_factor, :libflint), UInt64,
          (Ptr{nmod_poly_factor}, Ptr{nmod_poly}), &fac, &x)
  res = Array(Tuple{nmod_poly,Int}, fac._num)
  for i in 1:fac._num
    f = parent(x)()
    ccall((:nmod_poly_factor_get_nmod_poly, :libflint), Void,
            (Ptr{nmod_poly}, Ptr{nmod_poly_factor}, Int), &f, &fac, i-1)
    e = unsafe_load(fac.exp,i)
    res[i] = (f,e)
  end
  return res 
end  

function factor_squarefree(x::nmod_poly)
  fac = nmod_poly_factor(x._mod_n)
  ccall((:nmod_poly_factor_squarefree, :libflint), UInt64,
          (Ptr{nmod_poly_factor}, Ptr{nmod_poly}), &fac, &x)
  res = Array(Tuple{nmod_poly,Int}, fac._num)
  for i in 1:fac._num
    f = parent(x)()
    ccall((:nmod_poly_factor_get_nmod_poly, :libflint), Void,
            (Ptr{nmod_poly}, Ptr{nmod_poly_factor}, Int), &f, &fac, i-1)
    e = unsafe_load(fac.exp,i)
    res[i] = (f,e)
  end
  return res 
end  

function factor_distinct_deg(x::nmod_poly)
  !issquarefree(x) && error("Polynomial must be squarefree")
  degs = Array(Int, degree(x))
  degss = [ pointer(degs) ]
  fac = nmod_poly_factor(x._mod_n)
  ccall((:nmod_poly_factor_distinct_deg, :libflint), UInt64,
          (Ptr{nmod_poly_factor}, Ptr{nmod_poly}, Ptr{Ptr{Int}}),
          &fac, &x, degss)
  res = Array(Tuple{nmod_poly,Int}, fac._num)
  for i in 1:fac._num
    f = parent(x)()
    ccall((:nmod_poly_factor_get_nmod_poly, :libflint), Void,
            (Ptr{nmod_poly}, Ptr{nmod_poly_factor}, Int), &f, &fac, i-1)
    res[i] = (f,degs[i])
  end
  return res 
end  

function factor_shape(x::nmod_poly)
  res = Array(Int, degree(x))
  res2 = Array(Int, degree(x))
  res3 = Array(Tuple{Int, Int}, degree(x))
  k = Int(1)
  fill!(res,0)
  square_fac = factor_squarefree(x)
  for (f,i) in square_fac
    discdeg = factor_distinct_deg(f)
    for (g,j) in discdeg
      num = div(degree(g),j)
      res[j] += num*i
      res2[k] = j
      k += 1
    end
  end
  resize!(res2, k-1)
  res2 = unique(res2)
  resize!(res3, length(res2))
  k = Int(1)
  for j in 1:length(res2)
    res3[j] = (res2[j], res[res2[j]])
  end
  return res3
end  

################################################################################
#
#  Unsafe functions
#
################################################################################

function setcoeff!(x::nmod_poly, n::Int, y::UInt)
  ccall((:nmod_poly_set_coeff_ui, :libflint), Void, 
                   (Ptr{nmod_poly}, Int, UInt), &x, n, y)
end

function setcoeff!(x::nmod_poly, n::Int, y::fmpz)
  t = ZZ()
  ccall((:fmpz_mod_ui, :libflint), Void,
        (Ptr{fmpz}, Ptr{fmpz}, UInt64), &t, &y, parent(x)._n)
  u = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
  setcoeff!(x, n, u)
end

function setcoeff!(x::nmod_poly, n::Int, y::Integer)
  setcoeff!(x, n, ZZ(y))
end
  
function setcoeff!(x::nmod_poly, n::Int, y::Residue{fmpz})
  base_ring(x) != parent(y) && error("Incompatible parent objects")
  setcoeff!(x, n, y.data)
end

function add!(z::nmod_poly, x::nmod_poly, y::nmod_poly)
  ccall((:nmod_poly_add, :libflint), Void, 
          (Ptr{nmod_poly}, Ptr{nmod_poly},  Ptr{nmod_poly}), &z, &x, &y)
end

function sub!(z::nmod_poly, x::nmod_poly, y::nmod_poly)
  ccall((:nmod_poly_sub, :libflint), Void, 
          (Ptr{nmod_poly}, Ptr{nmod_poly},  Ptr{nmod_poly}), &z, &x, &y)
end

function mul!(z::nmod_poly, x::nmod_poly, y::nmod_poly)
  ccall((:nmod_poly_mul, :libflint), Void, 
          (Ptr{nmod_poly}, Ptr{nmod_poly},  Ptr{nmod_poly}), &z, &x, &y)
end

function _factor(x::nmod_poly)
  fac = nmod_poly_factor(x._mod_n)
  ccall((:nmod_poly_factor, :libflint), UInt64,
          (Ptr{nmod_poly_factor}, Ptr{nmod_poly}), &fac, &x)
  res = Array(Tuple{nmod_poly,Int}, fac._num)
  for i in 1:fac._num
    f = nmod_poly(x._mod_n)
    ccall((:nmod_poly_factor_get_nmod_poly, :libflint), Void,
            (Ptr{nmod_poly}, Ptr{nmod_poly_factor}, Int), &f, &fac, i-1)
    e = unsafe_load(fac.exp,i)
    res[i] = (f,e)
  end
  return res 
end  

function _lift(x::nmod_poly)
  z = fmpz_poly()
  ccall((:fmpz_poly_set_nmod_poly, :libflint), Void,
          (Ptr{fmpz_poly}, Ptr{nmod_poly}), &z, &x)
  return z
end

################################################################################
#
#  Parent object call overloads
#
################################################################################

function Base.call(R::NmodPolyRing)
  z = nmod_poly(R._n)
  z.parent = R
  return z
end

function Base.call(R::NmodPolyRing, x::fmpz)
  z = x*one(R)
  return z
end

function Base.call(R::NmodPolyRing, x::Integer)
  z = x*one(R)
  return z
end

function Base.call(R::NmodPolyRing, x::Residue{fmpz})
  base_ring(R) != parent(x) && error("Wrong parents")
  z = x*one(R)
  return z
end

function Base.call{T <: Integer}(R::NmodPolyRing, arr::Array{T, 1})
  length(arr) = 0 && error("Array must not be empty")
  z = nmod_poly(R._n, arr)
  z.parent = R
  return z
end

function Base.call(R::NmodPolyRing, arr::Array{fmpz, 1})
  length(arr) = 0 && error("Array must not be empty")
  z = nmod_poly(R._n, arr)
  z.parent = R
  return z
end

function Base.call(R::NmodPolyRing, arr::Array{Residue{fmpz}, 1})
  length(arr) == 0 && error("Array must have length > 0")
  (base_ring(R) != parent(arr[1])) && error("Wrong parents")
  z = nmod_poly(R._n, arr)
  z.parent = R
  return z
end

function Base.call(R::NmodPolyRing, x::fmpz_poly)
  z = nmod_poly(R._n, x)
  z.parent = R
  return z
end

################################################################################
#
#  Polynomial ring constructor
#
################################################################################

function PolynomialRing(R::ResidueRing{fmpz}, s::String)
  try parent_obj = NmodPolyRing(R, symbol(s))
    return parent_obj, parent_obj([R(0), R(1)])
  catch
    error("Not implemented (yet)")
  end
end

################################################################################
#
#  Parent checking
#
################################################################################

function check_parent(x::nmod_poly, y::nmod_poly)
  parent(x) != parent(y) && error("Parents must coincide")
  nothing
end

################################################################################
#
#  Canonicalization
#
################################################################################

canonical_unit(a::nmod_poly) = canonical_unit(lead(a))
