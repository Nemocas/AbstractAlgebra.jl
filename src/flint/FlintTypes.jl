###############################################################################
#
#   FlintTypes.jl : Parent and object types for Flint
#
###############################################################################

###############################################################################
#
#   FlintIntegerRing / fmpz
#
###############################################################################

type FlintIntegerRing <: Ring{Flint}
end

const FlintZZ = FlintIntegerRing()

type fmpz <: IntegerRingElem
    d::Int

    function fmpz()
        z = new()
        ccall((:fmpz_init, :libflint), Void, (Ptr{fmpz},), &z)
        finalizer(z, _fmpz_clear_fn)
        return z
    end

    function fmpz(x::Int)
        z = new()
        ccall((:fmpz_init_set_si, :libflint), Void, (Ptr{fmpz}, Int), &z, x)
        finalizer(z, _fmpz_clear_fn)
        return z
    end

    function fmpz(x::BigInt)
        z = new()
        ccall((:fmpz_init, :libflint), Void, (Ptr{fmpz},), &z)
        ccall((:fmpz_set_mpz, :libflint), Void, (Ptr{fmpz}, Ptr{BigInt}), &z, &x)
        finalizer(z, _fmpz_clear_fn)
        return z
    end

    function fmpz(x::Float64)
        !isinteger(x) && throw(InexactError())
        z = new()
        ccall((:fmpz_set_d, :libflint), Void, (Ptr{fmpz}, Cdouble), &z, x)
        finalizer(z, _fmpz_clear_fn)
        return z
    end

    fmpz(x::fmpz) = x
end

function _fmpz_clear_fn(a::fmpz)
   ccall((:fmpz_clear, :libflint), Void, (Ptr{fmpz},), &a)
end

###############################################################################
#
#   FlintRationalField / fmpq
#
###############################################################################

type FlintRationalField <: Field{Flint}
end

const FlintQQ = FlintRationalField()

type fmpq <: FractionElem{fmpz}
   num::Int
   den::Int

   function fmpq()
      z = new()
      ccall((:fmpq_init, :libflint), Void, (Ptr{fmpq},), &z)
      finalizer(z, _fmpq_clear_fn)
      return z
   end

   function fmpq(a::fmpz, b::fmpz)
      z = new()
      ccall((:fmpq_init, :libflint), Void, (Ptr{fmpq},), &z)
      ccall((:fmpq_set_fmpz_frac, :libflint), Void,
            (Ptr{fmpq}, Ptr{fmpz}, Ptr{fmpz}), &z, &a, &b)
      finalizer(z, _fmpq_clear_fn)
      return z
   end

   function fmpq(a::fmpz)
      z = new()
      ccall((:fmpq_init, :libflint), Void, (Ptr{fmpq},), &z)
      b = fmpz(1)
      ccall((:fmpq_set_fmpz_frac, :libflint), Void,
            (Ptr{fmpq}, Ptr{fmpz}, Ptr{fmpz}), &z, &a, &b)
      finalizer(z, _fmpq_clear_fn)
      return z
   end

   function fmpq(a::Int, b::Int)
      z = new()
      ccall((:fmpq_init, :libflint), Void, (Ptr{fmpq},), &z)
      ccall((:fmpq_set_si, :libflint), Void,
            (Ptr{fmpq}, Int, Int), &z, a, b)
      finalizer(z, _fmpq_clear_fn)
      return z
   end

   function fmpq(a::Int)
      z = new()
      ccall((:fmpq_init, :libflint), Void, (Ptr{fmpq},), &z)
      ccall((:fmpq_set_si, :libflint), Void,
            (Ptr{fmpq}, Int, Int), &z, a, 1)
      finalizer(z, _fmpq_clear_fn)
      return z
   end
end

_fmpq_clear_fn(a::fmpq) = ccall((:fmpq_clear, :libflint), Void, (Ptr{fmpq},), &a)

###############################################################################
#
#   FmpzPolyRing / fmpz_poly
#
###############################################################################

const FmpzPolyID = ObjectIdDict()

type FmpzPolyRing <: Ring{Flint}
   base_ring::Ring
   S::Symbol

   function FmpzPolyRing(s::Symbol)
      return try
         FmpzPolyID[s]
      catch
         FmpzPolyID[s] = new(FlintZZ, s)
      end
   end
end

type fmpz_poly <: PolyElem{fmpz}
   coeffs::Ptr{Void}
   alloc::Int
   length::Int
   parent::FmpzPolyRing

   function fmpz_poly()
      z = new()
      ccall((:fmpz_poly_init, :libflint), Void, (Ptr{fmpz_poly},), &z)
      finalizer(z, _fmpz_poly_clear_fn)
      return z
   end

   function fmpz_poly(a::Array{fmpz, 1})
      z = new()
      ccall((:fmpz_poly_init2, :libflint), Void, 
            (Ptr{fmpz_poly}, Int), &z, length(a))
      for i = 1:length(a)
         ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Void, 
                     (Ptr{fmpz_poly}, Int, Ptr{fmpz}), &z, i - 1, &a[i])
      end
      finalizer(z, _fmpz_poly_clear_fn)
      return z
   end

   function fmpz_poly(a::Int)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Void, (Ptr{fmpz_poly},), &z)
      ccall((:fmpz_poly_set_si, :libflint), Void, (Ptr{fmpz_poly}, Int), &z, a)
      finalizer(z, _fmpz_poly_clear_fn)
      return z
   end

   function fmpz_poly(a::fmpz)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Void, (Ptr{fmpz_poly},), &z)
      ccall((:fmpz_poly_set_fmpz, :libflint), Void, 
            (Ptr{fmpz_poly}, Ptr{fmpz}), &z, &a)
      finalizer(z, _fmpz_poly_clear_fn)
      return z
   end

   function fmpz_poly(a::fmpz_poly)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Void, (Ptr{fmpz_poly},), &z)
      ccall((:fmpz_poly_set, :libflint), Void, 
            (Ptr{fmpz_poly}, Ptr{fmpz_poly}), &z, &a)
      finalizer(z, _fmpz_poly_clear_fn)
      return z
   end
end

function _fmpz_poly_clear_fn(a::fmpz_poly)
   ccall((:fmpz_poly_clear, :libflint), Void, (Ptr{fmpz_poly},), &a)
end

###############################################################################
#
#   FmpqPolyRing / fmpq_poly
#
###############################################################################

const FmpqPolyID = ObjectIdDict()

type FmpqPolyRing <: Ring{Flint}
   base_ring::Field
   S::Symbol

   function FmpqPolyRing(R::FlintRationalField, s::Symbol)
      return try
         FmpqPolyID[s]
      catch
         FmpqPolyID[s] = new(R, s)
      end
   end
end

type fmpq_poly <: PolyElem{fmpq}
   coeffs::Ptr{Int}
   den::Int 
   alloc::Int
   length::Int
   parent::FmpqPolyRing

   function fmpq_poly()
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::Array{fmpq, 1})
      z = new()
      ccall((:fmpq_poly_init2, :libflint), Void, 
            (Ptr{fmpq_poly}, Int), &z, length(a))
      for i = 1:length(a)
         ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, 
                     (Ptr{fmpq_poly}, Int, Ptr{fmpq}), &z, i - 1, &a[i])
      end
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::Int)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_si, :libflint), Void, (Ptr{fmpq_poly}, Int), &z, a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::fmpz)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_fmpz, :libflint), Void, 
            (Ptr{fmpq_poly}, Ptr{fmpz}), &z, &a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::fmpq)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_fmpq, :libflint), Void,
            (Ptr{fmpq_poly}, Ptr{fmpq}), &z, &a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::fmpz_poly)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_fmpz_poly, :libflint), Void, 
            (Ptr{fmpq_poly}, Ptr{fmpz_poly}), &z, &a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::fmpq_poly)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set, :libflint), Void, 
            (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end
end

function _fmpq_poly_clear_fn(a::fmpq_poly)
   ccall((:fmpq_poly_clear, :libflint), Void, (Ptr{fmpq_poly},), &a)
end

###############################################################################
#
#   NmodPolyRing / nmod_poly
#
###############################################################################

const NmodPolyRingID = ObjectIdDict()

type NmodPolyRing <: Ring{Flint}
  base_ring::ResidueRing
  S::Symbol
  _n::UInt

  function NmodPolyRing(R::ResidueRing{fmpz}, s::Symbol)
    m = UInt(modulus(R))
    return try
       NmodPolyRingID[m, s]
    catch
       NmodPolyRingID[m, s] = new(R, s, m)
    end
  end
end

type nmod_poly <: PolyElem{Residue{fmpz}}
   _coeffs::Ptr{Void}
   _alloc::Int
   _length::Int
   _mod_n::UInt
   _mod_ninv::UInt
   _mod_norm::UInt
   parent::NmodPolyRing

   function nmod_poly(n::UInt)
      z = new()
      ccall((:nmod_poly_init, :libflint), Void, (Ptr{nmod_poly}, UInt), &z, n)
      finalizer(z, _nmod_poly_clear_fn)
      return z
   end

   function nmod_poly(n::UInt, a::UInt)
      z = new()
      ccall((:nmod_poly_init, :libflint), Void, (Ptr{nmod_poly}, UInt), &z, n)
      ccall((:nmod_poly_set_coeff_ui, :libflint), Void,
              (Ptr{nmod_poly}, Int, UInt), &z, 0, a)
      finalizer(z, _nmod_poly_clear_fn)
      return z
   end

   function nmod_poly(n::UInt, a::Int)
      z = new()
      ccall((:nmod_poly_init, :libflint), Void, (Ptr{nmod_poly}, UInt), &z, n)
      ccall((:nmod_poly_set_coeff_ui, :libflint), Void,
              (Ptr{nmod_poly}, Int, UInt), &z, 0, mod(a, n))
      finalizer(z, _nmod_poly_clear_fn)
      return z
   end

   function nmod_poly(n::UInt, arr::Array{fmpz, 1})
      z = new()
      ccall((:nmod_poly_init2, :libflint), Void,
            (Ptr{nmod_poly}, UInt, Int), &z, n, length(arr))
      for i in 1:length(arr)
         tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &arr[i])
         ccall((:nmod_poly_set_coeff_ui, :libflint), Void,
              (Ptr{nmod_poly}, Int, UInt), &z, i - 1, tt)
      end
      finalizer(z, _nmod_poly_clear_fn)
      return z
   end

   function nmod_poly(n::UInt, arr::Array{UInt, 1})
      z = new()
      ccall((:nmod_poly_init2, :libflint), Void,
            (Ptr{nmod_poly}, UInt, Int), &z, n, length(arr))
      for i in 1:length(arr)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Void,
              (Ptr{nmod_poly}, Int, UInt), &z, i - 1, arr[i])
      end
      finalizer(z, _nmod_poly_clear_fn)
      return z
   end

   function nmod_poly(n::UInt, arr::Array{Residue{fmpz}, 1})
      z = new()
      ccall((:nmod_poly_init2, :libflint), Void,
            (Ptr{nmod_poly}, UInt, Int), &z, n, length(arr))
      for i in 1:length(arr)
         tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &(arr[i]).data)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Void,
              (Ptr{nmod_poly}, Int, UInt), &z, i-1, tt)
      end
      finalizer(z, _nmod_poly_clear_fn)
      return z
   end

   function nmod_poly(n::UInt, f::fmpz_poly)
      z = new()
      ccall((:nmod_poly_init2, :libflint), Void,
            (Ptr{nmod_poly}, UInt, Int), &z, n, degree(f))
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
  exp::Ptr{Int} 
  _num::Int
  _alloc::Int
  _n::UInt
    
  function nmod_poly_factor(n::UInt)
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

###############################################################################
#
#   FmpzModPolyRing / fmpz_mod_poly
#
###############################################################################

const FmpzModPolyRingID = ObjectIdDict()

type FmpzModPolyRing <: Ring{Flint}
  base_ring::ResidueRing
  S::Symbol
  _n::fmpz

  function FmpzModPolyRing(R::ResidueRing{fmpz}, s::Symbol)
    m = modulus(R)
    return try
       FmpzModPolyRingID[m, s]
    catch
       FmpzModPolyRingID[m, s] = new(R, s, m)
    end
  end
end

type fmpz_mod_poly <: PolyElem{Residue{fmpz}}
   _coeffs::Ptr{Void}
   _alloc::Int
   _length::Int
   _p::Int
   parent::FmpzModPolyRing

   function fmpz_mod_poly(n::fmpz)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Void,
            (Ptr{fmpz_mod_poly}, Ptr{fmpz}), &z, &n)
      finalizer(z, _fmpz_mod_poly_clear_fn)
      return z
   end

   function fmpz_mod_poly(n::fmpz, a::UInt)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Void, 
            (Ptr{fmpz_mod_poly}, Ptr{fmpz}), &z, &n)
      ccall((:fmod_poly_set_coeff_ui, :libflint), Void,
              (Ptr{fmpz_mod_poly}, Int, UInt), &z, 0, a)
      finalizer(z, _fmpz_mod_poly_clear_fn)
      return z
   end

   function fmpz_mod_poly(n::fmpz, a::fmpz)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Void, 
            (Ptr{fmpz_mod_poly}, Ptr{fmpz}), &z, &n)
      ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Void,
              (Ptr{fmpz_mod_poly}, Int, Ptr{fmpz}), &z, 0, &a)
      finalizer(z, _fmpz_mod_poly_clear_fn)
      return z
   end

   function fmpz_mod_poly(n::fmpz, arr::Array{fmpz, 1})
      length(arr) == 0  && error("Array must have length > 0")
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Void,
            (Ptr{fmpz_mod_poly}, Ptr{fmpz}, Int), &z, &n, length(arr))
      for i in 1:length(arr)
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Void,
              (Ptr{fmpz_mod_poly}, Int, Ptr{fmpz}), &z, i - 1, &arr[i])
      end
      finalizer(z, _fmpz_mod_poly_clear_fn)
      return z
   end

   function fmpz_mod_poly(n::fmpz, arr::Array{Residue{fmpz}, 1})
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Void,
            (Ptr{fmpz_mod_poly}, Ptr{fmpz}, Int), &z, &n, length(arr))
      for i in 1:length(arr)
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Void,
              (Ptr{fmpz_mod_poly}, Int, Ptr{fmpz}), &z, i - 1, &(arr[i].data))
      end
      finalizer(z, _fmpz_mod_poly_clear_fn)
      return z
   end

   function fmpz_mod_poly(n::fmpz, f::fmpz_poly)
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Void,
            (Ptr{fmpz_mod_poly}, Ptr{fmpz}, Int), &z, &n, degree(f))
      ccall((:fmpz_mod_poly_set_fmpz_poly, :libflint), Void,
            (Ptr{fmpz_mod_poly}, Ptr{fmpz_poly}), &z, &f)
      finalizer(z, _fmpz_mod_poly_clear_fn)
      return z
   end

   function fmpz_mod_poly(n::fmpz, f::fmpz_mod_poly)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Void, (Ptr{fmpz_mod_poly}, Ptr{fmpz}), &z, &n)
      ccall((:fmpz_mod_poly_set, :libflint), Void,
            (Ptr{fmpz_mod_poly}, Ptr{fmpz_mod_poly}), &z, &f)
      finalizer(z, _fmpz_mod_poly_clear_fn)
      return z
   end
end

function _fmpz_mod_poly_clear_fn(x::fmpz_mod_poly)
  ccall((:fmpz_mod_poly_clear, :libflint), Void, (Ptr{fmpz_mod_poly}, ), &x)
end

type fmpz_mod_poly_factor
  poly::Ptr{fmpz_mod_poly}
  exp::Ptr{Int} 
  _num::Int
  _alloc::Int
  _n::fmpz
    
  function fmpz_mod_poly_factor(n::fmpz)
    z = new()
    ccall((:fmpz_mod_poly_factor_init, :libflint), Void,
            (Ptr{fmpz_mod_poly_factor}, ), &z)
    z._n = n
    finalizer(z, _fmpz_mod_poly_factor_clear_fn)
    return z
  end
end

function _fmpz_mod_poly_factor_clear_fn(a::fmpz_mod_poly_factor)
  ccall((:fmpz_mod_poly_factor_clear, :libflint), Void,
          (Ptr{fmpz_mod_poly_factor}, ), &a)
end

###############################################################################
#
#   FqNmodFiniteField / fq_nmod
#
###############################################################################

const FqNmodFiniteFieldID = Dict{Tuple{fmpz, Int, Symbol}, Field}()
const FqNmodFiniteFieldIDPol = Dict{Tuple{NmodPolyRing, nmod_poly, Symbol}, Field}()

type FqNmodFiniteField <: Field{Flint}
   p :: Int 
   n :: Int
   ninv :: Int
   norm :: Int
   sparse_modulus :: Int
   a :: Ptr{Void}
   j :: Ptr{Void}
   len :: Int
   mod_coeffs :: Ptr{Void}
   mod_alloc :: Int
   mod_length :: Int
   mod_n :: Int
   mod_ninv :: Int
   mod_norm :: Int
   inv_coeffs :: Ptr{Void}
   inv_alloc :: Int
   inv_length :: Int
   inv_n :: Int
   inv_ninv :: Int
   inv_norm :: Int
   var :: Ptr{Void}

   function FqNmodFiniteField(c::fmpz, deg::Int, s::Symbol)
      try
        return FqNmodFiniteFieldID[char, deg, s]
      catch
         d = new()
         ccall((:fq_nmod_ctx_init, :libflint), Void, 
               (Ptr{FqNmodFiniteField}, Ptr{fmpz}, Int, Ptr{UInt8}), 
			    &d, &c, deg, bytestring(string(s)))
         FqNmodFiniteFieldID[c, deg, s] = d
         finalizer(d, _FqNmodFiniteField_clear_fn)
         return d
      end
   end

   function FqNmodFiniteField(f::nmod_poly, s::Symbol)
      try
         return FqNmodFiniteFieldIDPol[parent(f), f, s]
      catch
         z = new()
         ccall((:fq_nmod_ctx_init_modulus, :libflint), Void, 
            (Ptr{FqNmodFiniteField}, Ptr{nmod_poly}, Ptr{UInt8}), 
	      &z, &f, bytestring(string(s)))
         FqNmodFiniteFieldIDPol[parent(f), f, s] = z
         finalizer(z, _FqNmodFiniteField_clear_fn)
         return z
      end
   end
end

function _FqNmodFiniteField_clear_fn(a :: FqNmodFiniteField)
   ccall((:fq_nmod_ctx_clear, :libflint), Void, (Ptr{FqNmodFiniteField},), &a)
end

type fq_nmod <: FiniteFieldElem
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   n :: Int
   ninv :: Int
   norm :: Int
   parent::FqNmodFiniteField

   function fq_nmod(ctx::FqNmodFiniteField)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      return d
   end

   function fq_nmod(ctx::FqNmodFiniteField, x::Int)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void,
            (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      ccall((:fq_nmod_set_si, :libflint), Void, 
                (Ptr{fq_nmod}, Int, Ptr{FqNmodFiniteField}), &d, x, &ctx)
      return d
   end

   function fq_nmod(ctx::FqNmodFiniteField, x::fmpz)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      ccall((:fq_nmod_set_fmpz, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{fmpz}, Ptr{FqNmodFiniteField}), &d, &x, &ctx)
      return d
   end

      function fq_nmod(ctx::FqNmodFiniteField, x::fq_nmod)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      ccall((:fq_nmod_set, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &x, &ctx)
      return d
   end
end

function _fq_nmod_clear_fn(a::fq_nmod)
   ccall((:fq_nmod_clear, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &a, &a.parent)
end

###############################################################################
#
#   FqFiniteField / fq
#
###############################################################################

const FqFiniteFieldID = Dict{Tuple{fmpz, Int, Symbol}, Field}()

const FqFiniteFieldIDPol = Dict{Tuple{fmpz_mod_poly, Symbol}, Field}()

type FqFiniteField <: Field{Flint}
   p::Int # fmpz
   sparse_modulus::Int
   a::Ptr{Void}
   j::Ptr{Void}
   len::Int
   mod_coeffs::Ptr{Void}
   mod_alloc::Int
   mod_length::Int
   mod_p::Int # fmpz
   inv_coeffs::Ptr{Void}
   inv_alloc::Int
   inv_length::Int
   inv_p::Int # fmpz
   var::Ptr{Void}

   function FqFiniteField(char::fmpz, deg::Int, s::Symbol)
      try
         return FqFiniteFieldID[char, deg, s]
      catch
         d = FqFiniteFieldID[char, deg, s] = new()
         finalizer(d, _FqFiniteField_clear_fn)
         ccall((:fq_ctx_init, :libflint), Void,
               (Ptr{FqFiniteField}, Ptr{fmpz}, Int, Ptr{UInt8}),
                  &d, &char, deg, bytestring(string(s)))
         return d
      end
   end
   
   function FqFiniteField(f::fmpz_mod_poly, s::Symbol)
      try
         return FqFiniteFieldIDPol[f, s]
      catch
         z = new()
         ccall((:fq_ctx_init_modulus, :libflint), Void,
               (Ptr{FqFiniteField}, Ptr{fmpz_mod_poly}, Ptr{UInt8}),
                  &z, &f, bytestring(string(s)))
         FqFiniteFieldIDPol[f, s] = z
         finalizer(z, _FqFiniteField_clear_fn)
         return z
      end
   end
end

function _FqFiniteField_clear_fn(a :: FqFiniteField)
   ccall((:fq_ctx_clear, :libflint), Void, (Ptr{FqFiniteField},), &a)
end

type fq <: FiniteFieldElem
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   parent::FqFiniteField

   function fq(ctx::FqFiniteField)
      d = new()
      ccall((:fq_init2, :libflint), Void, 
            (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      d.parent = ctx
      return d
   end

   function fq(ctx::FqFiniteField, x::Int)
      d = new()
      ccall((:fq_init2, :libflint), Void, 
            (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      ccall((:fq_set_si, :libflint), Void, 
                (Ptr{fq}, Int, Ptr{FqFiniteField}), &d, x, &ctx)
      d.parent = ctx
      return d
   end

   function fq(ctx::FqFiniteField, x::fmpz)
      d = new()
      ccall((:fq_init2, :libflint), Void, 
            (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      ccall((:fq_set_fmpz, :libflint), Void, 
            (Ptr{fq}, Ptr{fmpz}, Ptr{FqFiniteField}), &d, &x, &ctx)
      d.parent = ctx
      return d
   end

   function fq(ctx::FqFiniteField, x::fq)
      d = new()
      ccall((:fq_init2, :libflint), Void, 
            (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      ccall((:fq_set, :libflint), Void, 
            (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &d, &x, &ctx)
      d.parent = ctx
      return d
   end
end

function _fq_clear_fn(a::fq)
   ccall((:fq_clear, :libflint), Void, 
         (Ptr{fq}, Ptr{FqFiniteField}), &a, &a.parent)
end

###############################################################################
#
#   FlintPadicField / padic
#
###############################################################################

const PadicBase = ObjectIdDict()

type FlintPadicField <: Field{Flint}
   p::Int 
   pinv::Float64
   pow::Ptr{Void}
   minpre::Int
   maxpre::Int
   mode::Int
   prec_max::Int

   function FlintPadicField(p::fmpz, prec::Int)
      !isprime(p) && error("Prime base required in FlintPadicField")
      d = new()
      ccall((:padic_ctx_init, :libflint), Void, 
           (Ptr{FlintPadicField}, Ptr{fmpz}, Int, Int, Cint), 
                                     &d, &p, 0, 0, 1)
      finalizer(d, _padic_ctx_clear_fn)
      d.prec_max = prec
      return d
   end
end

function _padic_ctx_clear_fn(a::FlintPadicField)
   ccall((:padic_ctx_clear, :libflint), Void, (Ptr{FlintPadicField},), &a)
end

type padic <: PadicFieldElem
   u :: Int
   v :: Int
   N :: Int
   parent::FlintPadicField

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

###############################################################################
#
#   FmpzSeriesRing / fmpz_series
#
###############################################################################

const FmpzSeriesID = ObjectIdDict()

type FmpzSeriesRing <: Ring{Flint}
   base_ring::Ring
   prec_max::Int
   S::Symbol

   function FmpzSeriesRing(prec::Int, s::Symbol)
      return try
         FmpzSeriesID[prec, s]
      catch
         FmpzSeriesID[prec, s] = new(FlintZZ, prec, s)
      end
   end
end

type fmpz_series <: SeriesElem{fmpz}
   coeffs::Ptr{Void}
   alloc::Int
   length::Int
   prec :: Int
   parent::FmpzSeriesRing

   function fmpz_series()
      z = new()
      ccall((:fmpz_poly_init, :libflint), Void, 
            (Ptr{fmpz_series},), &z)
      finalizer(z, _fmpz_series_clear_fn)
      return z
   end
   
   function fmpz_series(a::Array{fmpz, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fmpz_poly_init2, :libflint), Void, 
            (Ptr{fmpz_series}, Int), &z, len)
      for i = 1:len
         ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Void, 
                     (Ptr{fmpz_series}, Int, Ptr{fmpz}), &z, i - 1, &a[i])
      end
      z.prec = prec
      finalizer(z, _fmpz_series_clear_fn)
      return z
   end
   
   function fmpz_series(a::fmpz_series)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Void, (Ptr{fmpz_series},), &z)
      ccall((:fmpz_poly_set, :libflint), Void, 
            (Ptr{fmpz_series}, Ptr{fmpz_series}), &z, &a)
      finalizer(z, _fmpz_series_clear_fn)
      return z
   end
end

function _fmpz_series_clear_fn(a::fmpz_series)
   ccall((:fmpz_poly_clear, :libflint), Void, (Ptr{fmpz_series},), &a)
end

###############################################################################
#
#   FmpqSeriesRing / fmpq_series
#
###############################################################################

const FmpqSeriesID = ObjectIdDict()

type FmpqSeriesRing <: Ring{Flint}
   base_ring::Ring
   prec_max::Int
   S::Symbol

   function FmpqSeriesRing(prec::Int, s::Symbol)
      return try
         FmpqSeriesID[prec, s]
      catch
         FmpqSeriesID[prec, s] = new(FlintQQ, prec, s)
      end
   end
end

type fmpq_series <: SeriesElem{fmpq}
   coeffs::Ptr{Void}
   den::Int
   alloc::Int
   length::Int
   prec :: Int
   parent::FmpqSeriesRing

   function fmpq_series()
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, 
            (Ptr{fmpq_series},), &z)
      finalizer(z, _fmpq_series_clear_fn)
      return z
   end
   
   function fmpq_series(a::Array{fmpq, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fmpq_poly_init2, :libflint), Void, 
            (Ptr{fmpq_series}, Int), &z, len)
      for i = 1:len
         ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, 
                     (Ptr{fmpq_series}, Int, Ptr{fmpq}), &z, i - 1, &a[i])
      end
      z.prec = prec
      finalizer(z, _fmpq_series_clear_fn)
      return z
   end
   
   function fmpq_series(a::fmpq_series)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_series},), &z)
      ccall((:fmpq_poly_set, :libflint), Void, 
            (Ptr{fmpq_series}, Ptr{fmpq_series}), &z, &a)
      finalizer(z, _fmpq_series_clear_fn)
      return z
   end
end

function _fmpq_series_clear_fn(a::fmpq_series)
   ccall((:fmpq_poly_clear, :libflint), Void, (Ptr{fmpq_series},), &a)
end

###############################################################################
#
#   FmpzModSeriesRing / fmpz_mod_series
#
###############################################################################

const FmpzModSeriesID = ObjectIdDict()

type FmpzModSeriesRing <: Ring{Flint}
   base_ring::ResidueRing
   prec_max::Int
   S::Symbol

   function FmpzModSeriesRing(R::Ring, prec::Int, s::Symbol)
      return try
         FmpzModSeriesID[R, prec, s]
      catch
         FmpzModSeriesID[R, prec, s] = new(R, prec, s)
      end
   end
end

type fmpz_mod_series <: SeriesElem{Residue{fmpz}}
   coeffs::Ptr{Void}
   alloc::Int
   length::Int
   p::Int
   prec :: Int
   parent::FmpzModSeriesRing

   function fmpz_mod_series(p::fmpz)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Void, 
            (Ptr{fmpz_mod_series}, Ptr{fmpz}), &z, &p)
      finalizer(z, _fmpz_mod_series_clear_fn)
      return z
   end
   
   function fmpz_mod_series(p::fmpz, a::Array{fmpz, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Void, 
            (Ptr{fmpz_mod_series}, Ptr{fmpz}, Int), &z, &p, len)
      for i = 1:len
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Void, 
                     (Ptr{fmpz_mod_series}, Int, Ptr{fmpz}), &z, i - 1, &a[i])
      end
      z.prec = prec
      finalizer(z, _fmpz_mod_series_clear_fn)
      return z
   end
   
   function fmpz_mod_series(p::fmpz, a::Array{Residue{fmpz}, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Void, 
            (Ptr{fmpz_mod_series}, Ptr{fmpz}, Int), &z, &p, len)
      for i = 1:len
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Void, 
                     (Ptr{fmpz_mod_series}, Int, Ptr{fmpz}), &z, i - 1, &data(a[i]))
      end
      z.prec = prec
      finalizer(z, _fmpz_mod_series_clear_fn)
      return z
   end
   
   function fmpz_mod_series(a::fmpz_mod_series)
      z = new()
      p = modulus(base_ring(parent(a)))
      ccall((:fmpz_mod_poly_init, :libflint), Void, 
            (Ptr{fmpz_mod_series}, Ptr{fmpz}), &z, &p)
      ccall((:fmpz_mod_poly_set, :libflint), Void, 
            (Ptr{fmpz_mod_series}, Ptr{fmpz_mod_series}), &z, &a)
      finalizer(z, _fmpz_mod_series_clear_fn)
      return z
   end
end

function _fmpz_mod_series_clear_fn(a::fmpz_mod_series)
   ccall((:fmpz_mod_poly_clear, :libflint), Void, (Ptr{fmpz_mod_series},), &a)
end

###############################################################################
#
#   FqSeriesRing / fq_series
#
###############################################################################

const FqSeriesID = ObjectIdDict()

type FqSeriesRing <: Ring{Flint}
   base_ring::Field
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

type fq_series <: SeriesElem{fq}
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

###############################################################################
#
#   FqNmodSeriesRing / fq_nmod_series
#
###############################################################################

const FqNmodSeriesID = ObjectIdDict()

type FqNmodSeriesRing <: Ring{Flint}
   base_ring::Field
   prec_max::Int
   S::Symbol

   function FqNmodSeriesRing(R::Ring{Flint}, prec::Int, s::Symbol)
      return try
         FqNmodSeriesID[R, prec, s]
      catch
         FqNmodSeriesID[R, prec, s] = new(R, prec, s)
      end
   end
end

type fq_nmod_series <: SeriesElem{fq_nmod}
   coeffs::Ptr{Void}
   alloc::Int
   length::Int
   prec :: Int
   parent::FqNmodSeriesRing

   function fq_nmod_series(ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_poly_init, :libflint), Void, 
            (Ptr{fq_nmod_series}, Ptr{FqNmodFiniteField}), &z, &ctx)
      finalizer(z, _fq_nmod_series_clear_fn)
      return z
   end
   
   function fq_nmod_series(ctx::FqNmodFiniteField, a::Array{fq_nmod, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fq_nmod_poly_init2, :libflint), Void, 
            (Ptr{fq_nmod_series}, Int, Ptr{FqNmodFiniteField}), &z, len, &ctx)
      for i = 1:len
         ccall((:fq_nmod_poly_set_coeff, :libflint), Void, 
               (Ptr{fq_nmod_series}, Int, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}),
                                               &z, i - 1, &a[i], &ctx)
      end
      z.prec = prec
      finalizer(z, _fq_nmod_series_clear_fn)
      return z
   end
   
   function fq_nmod_series(ctx::FqNmodFiniteField, a::fq_nmod_series)
      z = new()
      ccall((:fq_nmod_poly_init, :libflint), Void, 
            (Ptr{fq_nmod_series}, Ptr{FqNmodFiniteField}), &z, &ctx)
      ccall((:fq_nmod_poly_set, :libflint), Void, 
            (Ptr{fq_nmod_series}, Ptr{fq_nmod_series}, Ptr{FqNmodFiniteField}), &z, &a, &ctx)
      finalizer(z, _fq_nmod_series_clear_fn)
      return z
   end
end

function _fq_nmod_series_clear_fn(a::fq_nmod_series)
   ctx = base_ring(a)
   ccall((:fq_nmod_poly_clear, :libflint), Void,
         (Ptr{fq_nmod_series}, Ptr{FqNmodFiniteField}), &a, &ctx)
end

###############################################################################
#
#   FmpzMatSpace / fmpz_mat
#
###############################################################################

const FmpzMatID = ObjectIdDict()

# not really a mathematical ring
type FmpzMatSpace <: Ring{Flint}
   rows::Int
   cols::Int
   base_ring::FlintIntegerRing

   function FmpzMatSpace(r::Int, c::Int)
      return try
         FmpzMatID[r, c]
      catch
         FmpzMatID[r, c] = new(r, c, FlintZZ)
      end
   end
end

type fmpz_mat <: MatElem{fmpz}
   entries::Ptr{Void}
   r::Int
   c::Int
   rows::Ptr{Void}
   parent::Ring{Flint}

   # used by windows, not finalised!!
   function fmpz_mat()
      return new() 
   end

   function fmpz_mat(r::Int, c::Int)
      z = new()
      ccall((:fmpz_mat_init, :libflint), Void, 
            (Ptr{fmpz_mat}, Int, Int), &z, r, c)
      finalizer(z, _fmpz_mat_clear_fn)
      return z
   end

   function fmpz_mat(r::Int, c::Int, arr::Array{fmpz, 2})
      z = new()
      ccall((:fmpz_mat_init, :libflint), Void, 
            (Ptr{fmpz_mat}, Int, Int), &z, r, c)
      finalizer(z, _fmpz_mat_clear_fn)
      for i = 1:r
         for j = 1:c
            el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                       (Ptr{fmpz_mat}, Int, Int), &z, i - 1, j - 1)
            ccall((:fmpz_set, :libflint), Void,
                  (Ptr{fmpz}, Ptr{fmpz}), el, &arr[i, j])
         end
      end
      return z
   end

   function fmpz_mat{T <: Integer}(r::Int, c::Int, arr::Array{T, 2})
      z = new()
      ccall((:fmpz_mat_init, :libflint), Void, 
            (Ptr{fmpz_mat}, Int, Int), &z, r, c)
      finalizer(z, _fmpz_mat_clear_fn)
      for i = 1:r
         for j = 1:c
            el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                       (Ptr{fmpz_mat}, Int, Int), &z, i - 1, j - 1)
            ccall((:fmpz_set, :libflint), Void,
                  (Ptr{fmpz}, Ptr{fmpz}), el, &fmpz(arr[i, j]))
         end
      end
      return z
   end

   function fmpz_mat(r::Int, c::Int, d::fmpz)
      z = new()
      ccall((:fmpz_mat_init, :libflint), Void, 
            (Ptr{fmpz_mat}, Int, Int), &z, r, c)
      finalizer(z, _fmpz_mat_clear_fn)
      for i = 1:min(r, c)
         el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                    (Ptr{fmpz_mat}, Int, Int), &z, i - 1, i- 1)
         ccall((:fmpz_set, :libflint), Void,
               (Ptr{fmpz}, Ptr{fmpz}), el, &d)
      end
      return z
   end

   function fmpz_mat(m::fmpz_mat)
      z = new()
      ccall((:fmpz_mat_init_set, :libflint), Void, 
            (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &m)
      finalizer(z, _fmpz_mat_clear_fn)
      return z
   end
end

function _fmpz_mat_clear_fn(a::fmpz_mat)
   ccall((:fmpz_mat_clear, :libflint), Void, (Ptr{fmpz_mat},), &a)
end

###############################################################################
#
#   NmodMatSpace / nmod_mat
#
###############################################################################

const NmodMatID = ObjectIdDict()

type NmodMatSpace <: Ring{Flint}
  base_ring::ResidueRing
  _n::UInt
  rows::Int
  cols::Int

  function NmodMatSpace(R::ResidueRing{fmpz}, r::Int, c::Int)
    (r <= 0 || c <= 0) && error("Dimensions must be positive")
    fmpz(typemax(UInt)) < abs(R.modulus) &&
      error("Modulus of ResidueRing must less then ", fmpz(typemax(UInt)))
    try
      return NmodMatID[R, r, c]
    catch
      NmodMatID[R, r, c] = new(R, UInt(R.modulus), r, c)
    end
  end
end

type nmod_mat <: MatElem{Residue{fmpz}}
  entries::Ptr{Void}
  r::Int                  # Int
  c::Int                  # Int
  rows::Ptr{Void}
  _n::UInt                # mp_limb_t / Culong
  _ninv::UInt             # mp_limb_t / Culong
  _norm::UInt             # mp_limb_t / Culong
  parent::NmodMatSpace

  function nmod_mat(r::Int, c::Int, n::UInt)
    (r <= 0 || c <= 0) && error("Dimensions must be positive")
    z = new()
    ccall((:nmod_mat_init, :libflint), Void,
            (Ptr{nmod_mat}, Int, Int, UInt), &z, r, c, n)
    finalizer(z, _nmod_mat_clear_fn)
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{UInt, 2})
    (r <= 0 || c <= 0) && error("Dimensions must be positive")
    (size(arr) != (r,c)) && error("Array of wrong dimension")
    z = new()
    ccall((:nmod_mat_init, :libflint), Void,
            (Ptr{nmod_mat}, Int, Int, UInt), &z, r, c, n)
    finalizer(z, _nmod_mat_clear_fn)
    for i = 1:r
      for j = 1:c
        set_entry!(z, i, j, arr[i,j])
      end
    end
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{fmpz, 2})
    (r <= 0 || c <= 0) && error("Dimensions must be positive")
    (size(arr) != (r,c)) && error("Array of wrong dimension")
    z = new()
    ccall((:nmod_mat_init, :libflint), Void,
            (Ptr{nmod_mat}, Int, Int, UInt), &z, r, c, n)
    finalizer(z, _nmod_mat_clear_fn)
    for i = 1:r
      for j = 1:c
        set_entry!(z, i, j, arr[i,j])
      end
    end
    return z
  end

  function nmod_mat{T <: Integer}(r::Int, c::Int, n::UInt, arr::Array{T, 2})
    arr = map(fmpz, arr)
    return nmod_mat(r, c, n, arr)
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{Residue{fmpz}, 2})
    (r <= 0 || c <= 0) && error("Dimensions must be positive")
    (size(arr) != (r,c)) && error("Array of wrong dimension")
    z = new()
    ccall((:nmod_mat_init, :libflint), Void,
            (Ptr{nmod_mat}, Int, Int, UInt), &z, r, c, n)
    finalizer(z, _nmod_mat_clear_fn)
    for i = 1:r
      for j = 1:c
        set_entry!(z, i, j, arr[i,j])
      end
    end
    return z
  end

  function nmod_mat(n::UInt, b::fmpz_mat)
    z = new()
    ccall((:nmod_mat_init, :libflint), Void,
            (Ptr{nmod_mat}, Int, Int, UInt), &z, b.r, b.c, n)
    finalizer(z, _nmod_mat_clear_fn)
    ccall((:fmpz_mat_get_nmod_mat, :libflint), Void,
            (Ptr{nmod_mat}, Ptr{fmpz_mat}), &z, &b)
    return z
  end

  function nmod_mat(n::Int, b::fmpz_mat)
    (n < 0) && error("Modulus must be postive")
    return nmod_mat(UInt(n), b)
  end

  function nmod_mat(n::fmpz, b::fmpz_mat)
    (n < 0) && error("Modulus must be postive")
    (n > fmpz(typemax(UInt))) &&
          error("Exponent must be smaller then ", fmpz(typemax(UInt)))
    return nmod_mat(UInt(n), b) 
  end
end

function _nmod_mat_clear_fn(mat::nmod_mat)
  ccall((:nmod_mat_clear, :libflint), Void, (Ptr{nmod_mat}, ), &mat)
end

###############################################################################
#
#   FqPolyRing / fq_poly
#
###############################################################################

const FqPolyID = ObjectIdDict()

type FqPolyRing <: Ring{Flint}
   base_ring::Field
   S::Symbol

   function FqPolyRing(R::FqFiniteField, s::Symbol)
      if haskey(FqPolyID, (R, s))
         return FqPolyID[(R,s)]
      else
         z = new(R,s)
         FqPolyID[(R,s)] = z
         return z
      end
   end
end

type fq_poly <: PolyElem{fq}
   coeffs::Ptr{Void}
   alloc::Int
   length::Int
   parent::FqPolyRing

   function fq_poly()
      z = new()
      ccall((:fq_poly_init, :libflint), Void, (Ptr{fq_poly},), &z)   
      finalizer(z, _fq_poly_clear_fn)
      return z
   end

   function fq_poly(a::fq_poly)
      z = new()
      ctx = base_ring(parent(a))
      ccall((:fq_poly_init, :libflint), Void,
            (Ptr{fq_poly}, Ptr{FqFiniteField}), &z, &ctx)
      ccall((:fq_poly_set, :libflint), Void,
            (Ptr{fq_poly}, Ptr{fq_poly}, Ptr{FqFiniteField}),
            &z, &a, &ctx)
      finalizer(z, _fq_poly_clear_fn)
      return z
   end

   function fq_poly(a::fq)
      z = new()
      ctx = parent(a)
      ccall((:fq_poly_init, :libflint), Void,
            (Ptr{fq_poly}, Ptr{FqFiniteField}), &z, &ctx)
      ccall((:fq_poly_set_fq, :libflint), Void,
            (Ptr{fq_poly}, Ptr{fq}, Ptr{FqFiniteField}),
            &z, &a, &ctx)
      finalizer(z, _fq_poly_clear_fn)
      return z
   end

   function fq_poly(a::Array{fq, 1})
      z = new()
      ctx = parent(a[1])
      ccall((:fq_poly_init2, :libflint), Void,
            (Ptr{fq_poly}, Int, Ptr{FqFiniteField}),
            &z, length(a), &ctx)
      for i = 1:length(a)
         ccall((:fq_poly_set_coeff, :libflint), Void, 
               (Ptr{fq_poly}, Int, Ptr{fq}, Ptr{FqFiniteField}),
               &z, i - 1, &a[i], &ctx)
      end
      finalizer(z, _fq_poly_clear_fn)
      return z
   end

   function fq_poly(a::Array{fmpz, 1}, ctx::FqFiniteField)
      z = new()
      temp = ctx()
      ccall((:fq_poly_init2, :libflint), Void,
            (Ptr{fq_poly}, Int, Ptr{FqFiniteField}),
            &z, length(a), &ctx)
      for i = 1:length(a)
         temp = ctx(a[i])
         ccall((:fq_poly_set_coeff, :libflint), Void,
               (Ptr{fq_poly}, Int, Ptr{fq}, Ptr{FqFiniteField}),
               &z, i - 1, &temp, &ctx)
      end
      finalizer(z, _fq_poly_clear_fn)
      return z
   end

   function fq_poly(a::fmpz_poly, ctx::FqFiniteField)
      z = new()
      ccall((:fq_poly_init2, :libflint), Void,
            (Ptr{fq_poly}, Int, Ptr{FqFiniteField}),
            &z, length(a), &ctx)
      for i = 1:length(a)
         temp = ctx(coeff(a, i-1))
         ccall((:fq_poly_set_coeff, :libflint), Void,
               (Ptr{fq_poly}, Int, Ptr{fq}, Ptr{FqFiniteField}),
               &z, i - 1, &temp, &ctx)
      end
      finalizer(z, _fq_poly_clear_fn)
      return z
   end
end

function _fq_poly_clear_fn(a::fq_poly)
   ccall((:fq_poly_clear, :libflint), Void, (Ptr{fq_poly},), &a)
end

type fq_poly_factor
  poly::Ptr{fq_poly}
  exp::Ptr{Int} 
  num::Int
  alloc::Int
  base_field::FqFiniteField
    
  function fq_poly_factor(ctx::FqFiniteField)
    z = new()
    ccall((:fq_poly_factor_init, :libflint), Void,
         (Ptr{fq_poly_factor}, Ptr{FqFiniteField}), &z, &ctx)
    z.base_field = ctx
    finalizer(z, _fq_poly_factor_clear_fn)
    return z
  end
end

function _fq_poly_factor_clear_fn(a::fq_poly_factor)
   ccall((:fq_poly_factor_clear, :libflint), Void,
         (Ptr{fq_poly_factor}, Ptr{FqFiniteField}),
         &a, &(a.base_field))
end

###############################################################################
#
#   FqNmodPolyRing / fq_nmod_poly
#
###############################################################################

const FqNmodPolyID = ObjectIdDict()

type FqNmodPolyRing <: Ring{Flint}
   base_ring::Field
   S::Symbol

   function FqNmodPolyRing(R::FqNmodFiniteField, s::Symbol)
      if haskey(FqNmodPolyID, (R, s))
         return FqNmodPolyID[(R,s)]
      else
         z = new(R,s)
         FqNmodPolyID[(R,s)] = z
         return z
      end
   end
end

type fq_nmod_poly <: PolyElem{fq_nmod}
   coeffs::Ptr{Void}
   alloc::Int
   length::Int
   parent::FqNmodPolyRing

   function fq_nmod_poly()
      z = new()
      ccall((:fq_nmod_poly_init, :libflint), Void, (Ptr{fq_nmod_poly},), &z)   
      finalizer(z, _fq_nmod_poly_clear_fn)
      return z
   end

   function fq_nmod_poly(a::fq_nmod_poly)
      z = new()
      ctx = base_ring(parent(a))
      ccall((:fq_nmod_poly_init, :libflint), Void,
            (Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}), &z, &ctx)
      ccall((:fq_nmod_poly_set, :libflint), Void,
            (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
            &z, &a, &ctx)
      finalizer(z, _fq_nmod_poly_clear_fn)
      return z
   end

   function fq_nmod_poly(a::fq_nmod)
      z = new()
      ctx = parent(a)
      ccall((:fq_nmod_poly_init, :libflint), Void,
            (Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}), &z, &ctx)
      ccall((:fq_nmod_poly_set_fq_nmod, :libflint), Void,
            (Ptr{fq_nmod_poly}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}),
            &z, &a, &ctx)
      finalizer(z, _fq_nmod_poly_clear_fn)
      return z
   end

   function fq_nmod_poly(a::Array{fq_nmod, 1})
      z = new()
      ctx = parent(a[1])
      ccall((:fq_nmod_poly_init2, :libflint), Void,
            (Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}),
            &z, length(a), &ctx)
      for i = 1:length(a)
         ccall((:fq_nmod_poly_set_coeff, :libflint), Void, 
               (Ptr{fq_nmod_poly}, Int, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}),
               &z, i - 1, &a[i], &ctx)
      end
      finalizer(z, _fq_nmod_poly_clear_fn)
      return z
   end

   function fq_nmod_poly(a::Array{fmpz, 1}, ctx::FqNmodFiniteField)
      z = new()
      temp = ctx()
      ccall((:fq_nmod_poly_init2, :libflint), Void,
            (Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}),
            &z, length(a), &ctx)
      for i = 1:length(a)
         temp = ctx(a[i])
         ccall((:fq_nmod_poly_set_coeff, :libflint), Void,
               (Ptr{fq_nmod_poly}, Int, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}),
               &z, i - 1, &temp, &ctx)
      end
      finalizer(z, _fq_nmod_poly_clear_fn)
      return z
   end

   function fq_nmod_poly(a::fmpz_poly, ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_poly_init2, :libflint), Void,
            (Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}),
            &z, length(a), &ctx)
      for i = 1:length(a)
         temp = ctx(coeff(a,i-1))
         ccall((:fq_nmod_poly_set_coeff, :libflint), Void,
               (Ptr{fq_nmod_poly}, Int, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}),
               &z, i - 1, &temp, &ctx)
      end
      finalizer(z, _fq_nmod_poly_clear_fn)
      return z
   end
end

function _fq_nmod_poly_clear_fn(a::fq_nmod_poly)
   ccall((:fq_nmod_poly_clear, :libflint), Void, (Ptr{fq_nmod_poly},), &a)
end

type fq_nmod_poly_factor
  poly::Ptr{fq_nmod_poly}
  exp::Ptr{Int} 
  num::Int
  alloc::Int
  base_field::FqNmodFiniteField
    
  function fq_nmod_poly_factor(ctx::FqNmodFiniteField)
    z = new()
    ccall((:fq_nmod_poly_factor_init, :libflint), Void,
         (Ptr{fq_nmod_poly_factor}, Ptr{FqNmodFiniteField}), &z, &ctx)
    z.base_field = ctx
    finalizer(z, _fq_nmod_poly_factor_clear_fn)
    return z
  end
end

function _fq_nmod_poly_factor_clear_fn(a::fq_nmod_poly_factor)
   ccall((:fq_nmod_poly_factor_clear, :libflint), Void,
         (Ptr{fq_nmod_poly_factor}, Ptr{FqNmodFiniteField}),
         &a, &(a.base_field))
end

###############################################################################
#
#   FlintPermGroup / perm
#
###############################################################################

const FlintPermID = ObjectIdDict()

type FlintPermGroup <: Group{Flint}
   n::Int

   function FlintPermGroup(n::Int)
      return try
         FlintPermID[n]
      catch
         FlintPermID[n] = new(n)
      end
   end
end

type perm <: PermElem
   d::Array{Int, 1}
   parent::Group

   function perm(n::Int)
      p = new(Array(Int, n))
      ccall((:_perm_set_one, :libflint), Void,
            (Ref{Int}, Int), p.d, length(p.d))
      return p
   end

   function perm(a::Array{Int, 1})
      n = length(a)
      d = Array(Int, n)
      for i = 1:n
         d[i] = a[i] - 1
      end
      return new(d)
   end
end

   
