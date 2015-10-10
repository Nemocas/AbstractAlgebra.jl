###############################################################################
#
#   ArbTypes.jl : Parent and object types for Arb
#
#   Copyright (C) 2015 Tommy Hofmann
#   Copyright (C) 2015 Fredrik Johansson
#
###############################################################################

export ArbField, arb, AcbField, acb

arb_check_prec(p::Int) = (p >= 2 && p < (typemax(Int) >> 4)) || throw(ArgumentError("invalid precision"))

################################################################################
#
#  Types and memory management for ArbField
#
################################################################################

const ArbFieldID = ObjectIdDict()

type ArbField <: Field
  prec::Int
  
  function ArbField(p::Int = 256)
    arb_check_prec(p)
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
  size::UInt # mp_size_t
  d1::Int # mantissa_struct
  d2::Int
end

type mag_struct
  exp::Int # fmpz
  man::UInt # mp_limb_t
end

type arb_struct
  mid_exp::Int # fmpz
  mid_size::UInt # mp_size_t
  mid_d1::Int # mantissa_struct
  mid_d2::Int
  rad_exp::Int # fmpz
  rad_man::UInt
end

type arb <: FieldElem
  mid_exp::Int # fmpz
  mid_size::UInt # mp_size_t
  mid_d1::Int # mantissa_struct
  mid_d2::Int
  rad_exp::Int # fmpz
  rad_man::UInt
  parent::ArbField

  function arb()
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    finalizer(z, _arb_clear_fn)
    return z
  end

  function arb(x::Union{Int, UInt, Float64, fmpz, fmpq, BigFloat, AbstractString, arb}, p::Int)
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    _arb_set(z, x, p)
    finalizer(z, _arb_clear_fn)
    return z
  end

  function arb(x::Union{Int, UInt, Float64, fmpz, BigFloat, arb})
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    _arb_set(z, x)
    finalizer(z, _arb_clear_fn)
    return z
  end

  function arb(mid::arb, rad::arb)
    z = new()
    ccall((:arb_init, :libarb), Void, (Ptr{arb}, ), &z)
    ccall((:arb_set, :libarb), Void, (Ptr{arb}, Ptr{arb}), &z, &mid)
    ccall((:arb_add_error, :libarb), Void, (Ptr{arb}, Ptr{arb}), &z, &rad)
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

function deepcopy(a::arb)
  b = parent(a)()
  ccall((:arb_set, :libarb), Void, (Ptr{arb}, Ptr{arb}), &b, &a)
  return b
end

################################################################################
#
#  Types and memory management for AcbField
#
################################################################################

const AcbFieldID = ObjectIdDict()

type AcbField <: Field
  prec::Int
  
  function AcbField(p::Int = 256)
    arb_check_prec(p)
    try
      return AcbFieldID[p]::AcbField
    catch
      AcbFieldID[p] = new(p)::AcbField
    end
  end
end

prec(x::AcbField) = x.prec

type acb <: FieldElem
  real_mid_exp::Int     # fmpz
  real_mid_size::UInt # mp_size_t
  real_mid_d1::Int    # mantissa_struct
  real_mid_d2::Int
  real_rad_exp::Int     # fmpz
  real_rad_man::UInt
  imag_mid_exp::Int     # fmpz
  imag_mid_size::UInt # mp_size_t
  imag_mid_d1::Int    # mantissa_struct
  imag_mid_d2::Int
  imag_rad_exp::Int     # fmpz
  imag_rad_man::UInt
  parent::AcbField

  function acb()
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    finalizer(z, _acb_clear_fn)
    return z
  end

  function acb(x::Int)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_si, :libarb), Void, (Ptr{acb}, Int), &z, x)
    finalizer(z, _acb_clear_fn)
    return z
  end
  
  function acb(i::UInt)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_ui, :libarb), Void, (Ptr{acb}, UInt), &z, i)
    finalizer(z, _acb_clear_fn)
    return z
  end

  function acb(x::fmpz)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_fmpz, :libarb), Void, (Ptr{acb}, Ptr{fmpz}), &z, &x)
    finalizer(z, _acb_clear_fn)
    return z
  end

  function acb(x::arb)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_arb, :libarb), Void, (Ptr{acb}, Ptr{arb}), &z, &x)
    finalizer(z, _acb_clear_fn)
    return z
  end

  function acb(x::acb, p::Int)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_round, :libarb), Void,
                (Ptr{acb}, Ptr{acb}, Int), &z, &x, p)
    finalizer(z, _acb_clear_fn)
    return z
  end

  function acb(x::fmpz, p::Int)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_round_fmpz, :libarb), Void,
                (Ptr{acb}, Ptr{fmpz}, Int), &z, &x, p)
    finalizer(z, _acb_clear_fn)
    return z
  end

  function acb(x::fmpq, p::Int)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_fmpq, :libarb), Void,
                (Ptr{acb}, Ptr{fmpq}, Int), &z, &x, p)
    finalizer(z, _acb_clear_fn)
    return z
  end

  function acb(x::arb, p::Int)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_round_arb, :libarb), Void,
                (Ptr{acb}, Ptr{arb}, Int), &z, &x, p)
    finalizer(z, _acb_clear_fn)
    return z
  end

  function acb(x::Int, y::Int, p::Int)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_si_si, :libarb), Void,
                (Ptr{acb}, Int, Int, Int), &z, x, y, p)
    ccall((:acb_set_round, :libarb), Void,
                (Ptr{acb}, Ptr{acb}, Int), &z, &z, p)
    finalizer(z, _acb_clear_fn)
    return z
  end

  function acb(x::arb, y::arb, p::Int)
    z = new()
    ccall((:acb_init, :libarb), Void, (Ptr{acb}, ), &z)
    ccall((:acb_set_arb_arb, :libarb), Void,
                (Ptr{acb}, Ptr{arb}, Ptr{arb}, Int), &z, &x, &y, p)
    ccall((:acb_set_round, :libarb), Void,
                (Ptr{acb}, Ptr{acb}, Int), &z, &z, p)
    finalizer(z, _acb_clear_fn)
    return z
  end
end

function _acb_clear_fn(x::acb)
  ccall((:acb_clear, :libarb), Void, (Ptr{acb}, ), &x)
end

parent(x::acb) = x.parent

function deepcopy(a::acb)
  b = parent(a)()
  ccall((:acb_set, :libarb), Void, (Ptr{acb}, Ptr{acb}), &b, &a)
  return b
end

