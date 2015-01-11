###########################################################################################
#
#   fq.jl : Flint finite fields
#
###########################################################################################

export fq, FqFiniteField

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

FqFiniteFieldID = Dict{(BigInt, Int, Symbol), Field}()

type FqFiniteField{S} <: Field
   p :: Int # fmpz
   sparse_modulus :: Int
   a :: Ptr{Void}
   j :: Ptr{Void}
   len :: Int
   mod_coeffs :: Ptr{Void}
   mod_alloc :: Int
   mod_length :: Int
   mod_p :: Int # fmpz
   inv_coeffs :: Ptr{Void}
   inv_alloc :: Int
   inv_length :: Int
   inv_p :: Int # fmpz
   var :: Ptr{Void}

   function FqFiniteField(char::BigInt, deg::Int)
      try
         return FqFiniteFieldID[char, deg, S]
      catch
         d = FqFiniteFieldID[char, deg, S] = new()
         finalizer(d, _fq_ctx_clear_fn)
         temp = fmpz_readonly(char)
         ccall((:fq_ctx_init, :libflint), Void, 
               (Ptr{FqFiniteField}, Ptr{fmpz}, Int, Ptr{Uint8}), &d, &temp, deg, bytestring(string(S)))
         return d
      end
   end
end

function _fq_ctx_clear_fn(a :: FqFiniteField)
   ccall((:fq_ctx_clear, :libflint), Void, (Ptr{FqFiniteField},), &a)
end

type fq{S} <: FiniteFieldElem
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   parent::FqFiniteField{S}

   function fq(ctx::FqFiniteField)
      d = new()
      ccall((:fq_init2, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      return d
   end

   function fq(ctx::FqFiniteField, x::Int)
      z = fq{S}()
      ccall((:fq_set_si, :libflint), Void, 
                (Ptr{fq}, Int, Ptr{FqFiniteField}), &z, x, &ctx)
      return z
   end
   function fq(ctx::FqFiniteField, x::BigInt)
      z = fq{S}()
      temp = fmpz_readonly(x)
      ccall((:fq_set_fmpz, :libflint), Void, 
                (Ptr{fq}, Ptr{fmpz}, Ptr{FqFiniteField}), &z, &temp, &ctx)
      return z
   end

   fq(a::fq{S}) = a
end

function _fq_clear_fn{S}(a::fq{S})
   ccall((:fq_clear, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &a, &a.parent.ctx)
end

###########################################################################################
#
#   FiniteField constructor
#
###########################################################################################

function FiniteField(char::BigInt, deg::Int, s::String)
   S = symbol(s)
   parent_obj = FqFiniteField{S}(char, deg)

   return parent_obj # , parent_obj([ZZ(0), ZZ(1)]) 
end

FiniteField(char::Integer, deg::Int, s::String) = FiniteField(BigInt(char), deg, s)

