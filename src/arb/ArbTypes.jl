###############################################################################
#
#   ArbTypes.jl : Parent and object types for Arb
#
###############################################################################

export ArbField, arb

export radius, midpoint, zeropminf, indeterminate, contains, contains_zero,
       contains_negative, contains_positive, contains_nonnegative,
       contains_nonpositive, isnonzero, isexact, isint, ispositive,
       isnonnegative, isnegative, isnonpositive, sin, cos, tan, add!, mul!,
       sub!, div!, strongequal, prec, overlaps

################################################################################
#
#  Types and memory management
#
################################################################################

const ArbFieldID = ObjectIdDict()

type ArbField <: Field
  prec::Int
  
  function ArbField(p::Int = 256)
    try
      return ArbFieldID[p]::ArbField
    catch
      ArbFieldID[p] = new(p)
      return ArbFieldID[p]::ArbField
    end
  end
end

prec(x::ArbField) = x.prec

# these may be used for shallow operations
type arf_struct
  exp::Int # fmpz
  size::UInt64 # mp_size_t
  d1::Int64 # mantissa_struct
  d2::Int64
end

type mag_struct
  exp::Int # fmpz
  man::UInt64 # mp_limb_t
end

type arb_struct
  mid_exp::Int # fmpz
  mid_size::UInt64 # mp_size_t
  mid_d1::Int64 # mantissa_struct
  mid_d2::Int64
  rad_exp::Int # fmpz
  rad_man::UInt64
end

type arb <: FieldElem
  mid_exp::Int # fmpz
  mid_size::UInt64 # mp_size_t
  mid_d1::Int64 # mantissa_struct
  mid_d2::Int64
  rad_exp::Int # fmpz
  rad_man::UInt64
  parent::ArbField

  function arb()
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    finalizer(z, _arb_clear_fn)
    return z
  end

  function arb(x::arb, p::Int)
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    ccall((:arb_set_round, :libarb), Void,
                (Ptr{arb}, Ptr{arb}, Int), &z, &x, p)
    finalizer(z, _arb_clear_fn)
    return z
  end

  function arb(s::String, p::Int)
    s = bytestring(s)
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    err = ccall((:arb_set_str, :libarb), Int32, (Ptr{arb}, Ptr{Uint8}, Int), &z, s, p)
    finalizer(z, _arb_clear_fn)
    err == 0 || error("Invalid real string: $(repr(s))")
    return z
  end

  function arb(x::Int)
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    ccall((:arb_set_si, :libarb), Void, (Ptr{arb}, Int), &z, x)
    finalizer(z, _arb_clear_fn)
    return z
  end
 
  function arb(i::UInt)
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    ccall((:arb_set_ui, :libarb), Void, (Ptr{arb}, UInt), &z, i)
    finalizer(z, _arb_clear_fn)
    return z
  end

  function arb(x::fmpz)
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    ccall((:arb_set_fmpz, :libarb), Void, (Ptr{arb}, Ptr{fmpz}), &z, &x)
    finalizer(z, _arb_clear_fn)
    return z
  end
 
  function arb(x::fmpz, p::Int)
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    ccall((:arb_set_round_fmpz, :libarb), Void,
                (Ptr{arb}, Ptr{fmpz}, Int), &z, &x, p)
    finalizer(z, _arb_clear_fn)
    return z
  end
 
  function arb(x::fmpq, p::Int)
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    ccall((:arb_set_fmpq, :libarb), Void,
                (Ptr{arb}, Ptr{fmpq}, Int), &z, &x, p)
    finalizer(z, _arb_clear_fn)
    return z
  end
 
  #function arb(x::arf)
  #  z = new()
  #  ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
  #  ccall((:arb_set_arf, :libarb), Void, (Ptr{arb}, Ptr{arf}), &z, &x)
  #  finalizer(z, _arb_clear_fn)
  #  return z
  #end
end

function _arb_clear_fn(x::arb)
  ccall((:arb_clear, :libarb), Void, (Ptr{arb}, ), &x)
end

elem_type(x::ArbField) = arb

parent(x::arb) = x.parent

################################################################################
#
#  Parent object overloading
#
################################################################################

function call(r::ArbField)
  z = arb()
  z.parent = r
  return z
end

function call(r::ArbField, x::Int)
  z = arb(fmpz(x), r.prec)
  z.parent = r
  return z
end

function call(r::ArbField, x::UInt)
  z = arb(fmpz(x), r.prec)
  z.parent = r
  return z
end

function call(r::ArbField, x::fmpz)
  z = arb(x, r.prec)
  z.parent = r
  return z
end

function call(r::ArbField, x::fmpq)
  z = arb(x, r.prec)
  z.parent = r
  return z
end
  
#function call(r::ArbField, x::arf)
#  z = arb(arb(x), r.prec)
#  z.parent = r
#  return z
#end

#function call(r::ArbField, x::Float64)
#  return r(arf(x))
#end

function call(r::ArbField, x::arb)
  z = arb(x, r.prec)
  z.parent = r
  return z
end

function call(r::ArbField, x::String)
  z = arb(x, r.prec)
  z.parent = r
  return z
end

function call(r::ArbField, x::MathConst)
  if x == pi
    z = pi_arb(r.prec)
    z.parent = r
    return z
  elseif x == e
    z = e_arb(r.prec)
    z.parent = r 
    return z
  else
    error("constant not supported")
  end
end

################################################################################
#
#  String I/O
#
################################################################################

function show(io::IO, x::arb)
  d = ceil(parent(x).prec * 0.30102999566398119521)
  cstr = ccall((:arb_get_str, :libarb), Ptr{UInt8}, (Ptr{arb}, Int, UInt),
                                                  &x, Int(d), UInt(0))
  print(io, bytestring(cstr))
  ccall((:flint_free, :libflint), Void, (Ptr{UInt8},), cstr)
end

################################################################################
#
#  Comparison and containment
#
################################################################################

function strongequal(x::arb, y::arb)
  r = ccall((:arb_equal, :libarb), Cint, (Ptr{arb}, Ptr{arb}), &x, &y)
  return Bool(r)
end

function overlaps(x::arb, y::arb)
  r = ccall((:arb_overlaps, :libarb), Cint, (Ptr{arb}, Ptr{arb}), &x, &y)
  return Bool(r)
end

#function contains(x::arb, y::arf)
#  r = ccall((:arb_contains_arf, :libarb), Cint, (Ptr{arb}, Ptr{arf}), &x, &y)
#  return Bool(r)
#end

function contains(x::arb, y::fmpq)
  r = ccall((:arb_contains_fmpq, :libarb), Cint, (Ptr{arb}, Ptr{fmpq}), &x, &y)
  return Bool(r)
end

function contains(x::arb, y::fmpz)
  r = ccall((:arb_contains_fmpz, :libarb), Cint, (Ptr{arb}, Ptr{fmpz}), &x, &y)
  return Bool(r)
end

function contains(x::arb, y::Int)
  r = ccall((:arb_contains_si, :libarb), Cint, (Ptr{arb}, Ptr{fmpz}), &x, y)
  return Bool(r)
end

function contains(x::arb, y::BigFloat)
  r = ccall((:arb_contains_mpfr, :libarb), Cint,
              (Ptr{arb}, Ptr{BigFloat}), &x, &y)
  return Bool(r)
end

function contains(x::arb, y::arb)
  r = ccall((:arb_contains, :libarb), Cint, (Ptr{arb}, Ptr{arb}), &x, &y)
  return Bool(r)
end

for (s,f) in (("contains_zero", "arb_contains_zero"),
              ("contains_negative", "arb_contains_negative"),
              ("contains_positive", "arb_contains_positive"),
              ("contains_nonpositive", "arb_contains_nonpositive"),
              ("contains_nonnegative", "arb_contains_nonnegative"))
  @eval begin
    function($(symbol(s)))(x::arb)
      r = ccall(($f, :libarb), Cint, (Ptr{arb}, ), &x)
      return Bool(r)
    end
  end
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::arb, y::arb)
    return Bool(ccall((:arb_eq, :libarb), Cint, (Ptr{arb}, Ptr{arb}), &x, &y))
end

function !=(x::arb, y::arb)
    return Bool(ccall((:arb_ne, :libarb), Cint, (Ptr{arb}, Ptr{arb}), &x, &y))
end

function >(x::arb, y::arb)
    return Bool(ccall((:arb_gt, :libarb), Cint, (Ptr{arb}, Ptr{arb}), &x, &y))
end

function >=(x::arb, y::arb)
    return Bool(ccall((:arb_ge, :libarb), Cint, (Ptr{arb}, Ptr{arb}), &x, &y))
end

function <(x::arb, y::arb)
    return Bool(ccall((:arb_lt, :libarb), Cint, (Ptr{arb}, Ptr{arb}), &x, &y))
end

function <=(x::arb, y::arb)
    return Bool(ccall((:arb_le, :libarb), Cint, (Ptr{arb}, Ptr{arb}), &x, &y))
end

==(x::arb, y::Int) = x == arb(y)
!=(x::arb, y::Int) = x != arb(y)
<=(x::arb, y::Int) = x <= arb(y)
>=(x::arb, y::Int) = x >= arb(y)
<(x::arb, y::Int) = x < arb(y)
>(x::arb, y::Int) = x > arb(y)

==(x::Int, y::arb) = arb(x) == y
!=(x::Int, y::arb) = arb(x) != y
<=(x::Int, y::arb) = arb(x) <= y
>=(x::Int, y::arb) = arb(x) >= y
<(x::Int, y::arb) = arb(x) < y
>(x::Int, y::arb) = arb(x) > y

==(x::arb, y::fmpz) = x == arb(y)
!=(x::arb, y::fmpz) = x != arb(y)
<=(x::arb, y::fmpz) = x <= arb(y)
>=(x::arb, y::fmpz) = x >= arb(y)
<(x::arb, y::fmpz) = x < arb(y)
>(x::arb, y::fmpz) = x > arb(y)

==(x::fmpz, y::arb) = arb(x) == y
!=(x::fmpz, y::arb) = arb(x) != y
<=(x::fmpz, y::arb) = arb(x) <= y
>=(x::fmpz, y::arb) = arb(x) >= y
<(x::fmpz, y::arb) = arb(x) < y
>(x::fmpz, y::arb) = arb(x) > y

################################################################################
#
#  Predicates
#
################################################################################

for (s,f) in (("iszero", "arb_is_zero"),
              ("isnonzero", "arb_is_nonzero"),
              ("isone", "arb_is_one"),
              ("isfinite", "arb_is_finite"),
              ("isexact", "arb_is_exact"),
              ("isint", "arb_is_int"),
              ("ispositive", "arb_is_positive"),
              ("isnonnegative", "arb_is_nonnegative"),
              ("isnegative", "arb_is_negative"),
              ("isnonpositive", "arb_is_nonpositive"))
  @eval begin
    function($(symbol(s)))(x::arb)
      i = ccall(($f, :libarb), Cint, (Ptr{arb},), &x)
      return Bool(i)
    end
  end
end

function radius(x::arb)
  t = mag_struct(x.rad_exp, x.rad_man)
  u = arf_struct(0, 0, 0, 0)
  z = parent(x)()
  ccall((:arf_set_mag, :libarb), Void, (Ptr{arf_struct}, Ptr{mag_struct}), &u, &t)
  z.mid_exp = u.exp
  z.mid_size = u.size
  z.mid_d1 = u.d1
  z.mid_d2 = u.d2
  return z
end

function midpoint(x::arb)
  t = arf_struct(x.mid_exp, x.mid_size, x.mid_d1, x.mid_d2)
  u = arf_struct(0, 0, 0, 0)
  z = parent(x)()
  ccall((:arf_set, :libarb), Void, (Ptr{arf_struct}, Ptr{arf_struct}), &u, &t)
  z.mid_exp = u.exp
  z.mid_size = u.size
  z.mid_d1 = u.d1
  z.mid_d2 = u.d2
  return z
end

for (s,f) in (("iszero", "arb_is_zero"), ("isnonzero", "arb_is_nonzero"),
              ("isone", "arb_is_one"), ("isfinite", "arb_is_finite"),
              ("isexact", "arb_is_exact"), ("isint", "arb_is_int"),
              ("ispositive", "arb_is_positive"),
              ("isnonnegative", "arb_is_nonnegative"),
              ("isnegative", "arb_is_negative"),
              ("isnonnegative", "arb_is_nonnegative"))
  @eval begin
    function($(symbol(s)))(x::arb)
      return Bool(ccall(($f, :libarb), Cint, (Ptr{arb}, ), &x))
    end
  end
end

################################################################################
#
#  Unary operations
#
################################################################################

function abs(x::arb)
  z = parent(x)()
  ccall((:arb_abs, :libarb), Void, (Ptr{arb}, Ptr{arb}), &z, &x)
  return z
end

################################################################################
#
#  Binary operations
#
################################################################################

function -(x::arb)
  z = parent(x)()
  ccall((:arb_neg_round, :libarb), Void, (Ptr{arb}, Ptr{arb}), &z, &x)
  return z
end

for (s,f) in ((:+,"arb_add"), (:*,"arb_mul"), (:/, "arb_div"), (:-,"arb_sub"))
  @eval begin
    function ($s)(x::arb, y::arb)
      z = parent(x)()
      ccall(($f, :libarb), Void, (Ptr{arb}, Ptr{arb}, Ptr{arb}, Int),
                           &z, &x, &y, parent(x).prec)
      return z
    end
  end
end

for (f,s) in ((:+, "add"), (:*, "mul"))
  @eval begin
    #function ($f)(x::arb, y::arf)
    #  z = parent(x)()
    #  ccall(($("arb_"*s*"_arf"), :libarb), Void,
    #              (Ptr{arb}, Ptr{arb}, Ptr{arf}, Int),
    #              &z, &x, &y, parent(x).prec)
    #  return z
    #end

    #($f)(x::arf, y::arb) = ($f)(y, x)

    function ($f)(x::arb, y::UInt)
      z = parent(x)()
      ccall(($("arb_"*s*"_ui"), :libarb), Void,
                  (Ptr{arb}, Ptr{arb}, UInt, Int),
                  &z, &x, y, parent(x).prec)
      return z
    end

    ($f)(x::UInt, y::arb) = ($f)(y, x)

    function ($f)(x::arb, y::Int)
      z = parent(x)()
      ccall(($("arb_"*s*"_si"), :libarb), Void,
      (Ptr{arb}, Ptr{arb}, Int, Int), &z, &x, y, parent(x).prec)
      return z
    end

    ($f)(x::Int, y::arb) = ($f)(y,x)

    function ($f)(x::arb, y::fmpz)
      z = parent(x)()
      ccall(($("arb_"*s*"_fmpz"), :libarb), Void,
                  (Ptr{arb}, Ptr{arb}, Ptr{fmpz}, Int),
                  &z, &x, &y, parent(x).prec)
      return z
    end

    ($f)(x::fmpz, y::arb) = ($f)(y,x)
  end
end

#function -(x::arb, y::arf)
#  z = parent(x)()
#  ccall((:arb_sub_arf, :libarb), Void,
#              (Ptr{arb}, Ptr{arb}, Ptr{arf}, Int), &z, &x, &y, parent(x).prec)
#  return z
#end

#-(x::arf, y::arb) = -(y - x)

function -(x::arb, y::UInt)
  z = parent(x)()
  ccall((:arb_sub_ui, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, UInt, Int), &z, &x, y, parent(x).prec)
  return z
end

-(x::UInt, y::arb) = -(y - x)

function -(x::arb, y::Int)
  z = parent(x)()
  ccall((:arb_sub_si, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Int, Int), &z, &x, y, parent(x).prec)
  return z
end

-(x::Int, y::arb) = -(y - x)

function -(x::arb, y::fmpz)
  z = parent(x)()
  ccall((:arb_sub_fmpz, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Ptr{fmpz}, Int),
              &z, &x, &y, parent(x).prec)
  return z
end

-(x::fmpz, y::arb) = -(y-x)

#function /(x::arb, y::arf)
#  z = parent(x)()
#  ccall((:arb_div_arf, :libarb), Void,
#              (Ptr{arb}, Ptr{arb}, Ptr{arf}, Int), &z, &x, &y, parent(x).prec)
#  return z
#end

function /(x::arb, y::UInt)
  z = parent(x)()
  ccall((:arb_div_ui, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, UInt, Int), &z, &x, y, parent(x).prec)
  return z
end

function /(x::arb, y::Int)
  z = parent(x)()
  ccall((:arb_div_si, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Int, Int), &z, &x, y, parent(x).prec)
  return z
end

function /(x::arb, y::fmpz)
  z = parent(x)()
  ccall((:arb_div_fmpz, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Ptr{fmpz}, Int),
              &z, &x, &y, parent(x).prec)
  return z
end

function ^(x::arb, y::arb)
  z = parent(x)()
  ccall((:arb_pow, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Ptr{arb}, Int), &z, &x, &y, parent(x).prec)
  return z
end

function ^(x::arb, y::fmpz)
  z = parent(x)()
  ccall((:arb_pow_fmpz, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Ptr{fmpz}, Int),
              &z, &x, &y, parent(x).prec)
  return z
end

^(x::arb, y::Int) = ^(x, fmpz(y))

function ^(x::arb, y::UInt)
  z = parent(x)()
  ccall((:arb_pow_ui, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, UInt, Int), &z, &x, y, parent(x).prec)
  return z
end

function inv(x::arb)
  z = parent(x)()
  ccall((:arb_inv, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Int), &z, &x, parent(x).prec)
  return parent(x)(z)
end

################################################################################
#
#  Unsafe binary operations
#
################################################################################

for (s,f) in (("add!","arb_add"), ("mul!","arb_mul"), ("div!", "arb_div"),
              ("sub!","arb_sub"))
  @eval begin
    function ($(symbol(s)))(z::arb, x::arb, y::arb)
      ccall(($f, :libarb), Void, (Ptr{arb}, Ptr{arb}, Ptr{arb}, Int),
                           &z, &x, &y, parent(x).prec)
    end
  end
end

################################################################################
#
#  Real valued functions
#
################################################################################

function log(x::arb)
  z = parent(x)()
  ccall((:arb_log, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Int), &z, &x, parent(x).prec)
  return z
end

function exp(x::arb)
  z = parent(x)()
  ccall((:arb_exp, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Int), &z, &x, parent(x).prec)
  return z
end

function sqrt(x::arb)
  z = parent(x)()
  ccall((:arb_sqrt, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Int), &z, &x, parent(x).prec)
  return z
end

function sin(x::arb)
  z = parent(x)()
  ccall((:arb_sin, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Int), &z, &x, parent(x).prec)
  return z
end

function cos(x::arb)
  z = parent(x)()
  ccall((:arb_cos, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Int), &z, &x, parent(x).prec)
  return z
end

function tan(x::arb)
  z = parent(x)()
  ccall((:arb_tan, :libarb), Void,
              (Ptr{arb}, Ptr{arb}, Int), &z, &x, parent(x).prec)
  return z
end

################################################################################
#
#  Constants
#
################################################################################

function pi_arb(p::Int)
  z = ArbField(p)()
  ccall((:arb_const_pi, :libarb), Void, (Ptr{arb}, Int), &z, p)
  return z
end

function e_arb(p::Int)
  z = ArbField(p)()
  ccall((:arb_const_e, :libarb), Void, (Ptr{arb}, Int), &z, p)
  return z
end

