###########################################################################################
#
#   fmpz_mat.jl : Flint matrices over ZZ
#
###########################################################################################

import Base: getindex, setindex!, transpose

export fmpz_mat, MatrixSpace, getindex, setindex!, rows, cols, charpoly, determinant,
       determinant_divisor, determinant_given_divisor, gram, hadamard, is_hadamard, hnf,
       is_hnf, hnf_with_transform, hnf_modular, lll, lll_with_transform, lll_with_removal,
       lll_with_removal_transform, nullspace, rank, rref, reduce_mod, snf, snf_diagonal,
       is_snf, solve, solve_dixon, trace, transpose

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

FmpzMatID = ObjectIdDict()

# not really a mathematical ring
type FmpzMatSpace <: Ring
   rows::Int
   cols::Int
   base_ring::IntegerRing

   function FmpzMatSpace(r::Int, c::Int)
      return try
         FmpzMatID[r, c]
      catch
         FmpzMatID[r, c] = new(r, c, ZZ)
      end
   end
end

type fmpz_mat <: MatElem
   entries::Ptr{Void}
   r::Int
   c::Int
   rows::Ptr{Void}
   parent

   function fmpz_mat(r::Int, c::Int)
      z = new()
      ccall((:fmpz_mat_init, :libflint), Void, (Ptr{fmpz_mat}, Int, Int), &z, r, c)
      finalizer(z, _fmpz_mat_clear_fn)
      return z
   end

   function fmpz_mat(r::Int, c::Int, arr::Array{BigInt, 2})
      z = new()
      ccall((:fmpz_mat_init, :libflint), Void, (Ptr{fmpz_mat}, Int, Int), &z, r, c)
      finalizer(z, _fmpz_mat_clear_fn)
      for i = 1:r
         for j = 1:c
            el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                       (Ptr{fmpz_mat}, Int, Int), &z, i - 1, j - 1)
            ccall((:fmpz_set_mpz, :libflint), Void,
                  (Ptr{fmpz}, Ptr{BigInt}), el, &arr[i, j])
         end
      end
      return z
   end

   function fmpz_mat{T <: Integer}(r::Int, c::Int, arr::Array{T, 2})
      z = new()
      ccall((:fmpz_mat_init, :libflint), Void, (Ptr{fmpz_mat}, Int, Int), &z, r, c)
      finalizer(z, _fmpz_mat_clear_fn)
      for i = 1:r
         for j = 1:c
            el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                       (Ptr{fmpz_mat}, Int, Int), &z, i - 1, j - 1)
            ccall((:fmpz_set_mpz, :libflint), Void,
                  (Ptr{fmpz}, Ptr{BigInt}), el, &BigInt(arr[i, j]))
         end
      end
      return z
   end

   function fmpz_mat(r::Int, c::Int, d::BigInt)
      z = new()
      ccall((:fmpz_mat_init, :libflint), Void, (Ptr{fmpz_mat}, Int, Int), &z, r, c)
      finalizer(z, _fmpz_mat_clear_fn)
      for i = 1:min(r, c)
         el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                    (Ptr{fmpz_mat}, Int, Int), &z, i - 1, i- 1)
         ccall((:fmpz_set_mpz, :libflint), Void,
               (Ptr{fmpz}, Ptr{BigInt}), el, &d)
      end
      return z
   end

   function fmpz_mat(m::fmpz_mat)
      z = new()
      ccall((:fmpz_mat_init_set, :libflint), Void, (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &m)
      finalizer(z, _fmpz_mat_clear_fn)
      return z
   end
end

function _fmpz_mat_clear_fn(a::fmpz_mat)
   ccall((:fmpz_mat_clear, :libflint), Void, (Ptr{fmpz_mat},), &a)
end

elem_type(::FmpzMatSpace) = fmpz_mat

base_ring(a::FmpzMatSpace) = a.base_ring

parent(a::fmpz_mat) = a.parent

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function getindex(a::fmpz_mat, r::Int, c::Int)
   (r > a.parent.rows || c > a.parent.cols) && throw(BoundsError())
   z = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
             (Ptr{fmpz_mat}, Int, Int), &a, r - 1, c - 1)
   u = BigInt()
   ccall((:fmpz_get_mpz, :libflint), Void, (Ptr{BigInt}, Ptr{fmpz}), &u, z)
   return u
end
 
function setindex!(a::fmpz_mat, d::BigInt, r::Int, c::Int)
   (r > a.parent.rows || c > a.parent.cols) && throw(BoundsError())
   z = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
             (Ptr{fmpz_mat}, Int, Int), &a, r - 1, c - 1)
   ccall((:fmpz_set_mpz, :libflint), Void, (Ptr{fmpz}, Ptr{BigInt}), z, &d)
end

setindex!(a::fmpz_mat, d::Integer, r::Int, c::Int) = setindex!(a, ZZ(d), r, c)

rows(a::fmpz_mat) = a.parent.rows

cols(a::fmpz_mat) = a.parent.cols

zero(a::FmpzMatSpace) = a()

one(a::FmpzMatSpace) = a(1)

iszero(a::fmpz_mat) = ccall((:fmpz_mat_is_zero, :libflint), Bool, (Ptr{fmpz_mat},), &a)

isone(a::fmpz_mat) = ccall((:fmpz_mat_is_one, :libflint), Bool, (Ptr{fmpz_mat},), &a)

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit(a::fmpz_mat) = canonical_unit(a[0, 0])

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, a::FmpzMatSpace)
   print(io, "Matrix Space of ", a.rows, " rows and ", a.cols, " columns over ")
   print(io, "Integer Ring")
end

function show(io::IO, a::fmpz_mat)
   rows = a.parent.rows
   cols = a.parent.cols
   for i = 1:rows
      print(io, "[")
      for j = 1:cols
         print(io, a[i, j])
         if j != cols
            print(io, " ")
         end
      end
      print(io, "]")
      if i != rows
         println("")
      end
   end
end

show_minus_one(::Type{fmpz_mat}) = show_minus_one(BigInt)

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -(x::fmpz_mat)
   z = parent(x)()
   ccall((:fmpz_mat_neg, :libflint), Void, (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &x)
   return z
end

function transpose(x::fmpz_mat)
   z = parent(x)()
   ccall((:fmpz_mat_transpose, :libflint), Void, (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &x)
   return z
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

function +(x::fmpz_mat, y::fmpz_mat)
   (rows(x) != rows(y) || cols(x) != cols(y)) && error("Incompatible matrix dimensions")
   z = parent(x)()
   ccall((:fmpz_mat_add, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat},  Ptr{fmpz_mat}), 
               &z, &x, &y)
   return z
end

function -(x::fmpz_mat, y::fmpz_mat)
   (rows(x) != rows(y) || cols(x) != cols(y)) && error("Incompatible matrix dimensions")
   z = parent(x)()
   ccall((:fmpz_mat_sub, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat},  Ptr{fmpz_mat}), 
               &z, &x, &y)
   return z
end

function *(x::fmpz_mat, y::fmpz_mat)
   cols(x) != rows(y) && error("Incompatible matrix dimensions")
   if rows(x) == cols(y) && rows(x) == cols(x)
      parz = parent(x)
   else
      parz = FmpzMatSpace(rows(x), cols(y))
   end
   z = parz()
   ccall((:fmpz_mat_mul, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat},  Ptr{fmpz_mat}), 
               &z, &x, &y)
   return z
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *(x::Int, y::fmpz_mat)
   z = parent(y)()
   ccall((:fmpz_mat_scalar_mul_si, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Int), &z, &y, x)
   return z
end

function *(x::BigInt, y::fmpz_mat)
   z = parent(y)()
   temp = fmpz_readonly(x)
   ccall((:fmpz_mat_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{fmpz_readonly}), &z, &y, &temp)
   return z
end

*(x::fmpz_mat, y::Int) = y*x

*(x::fmpz_mat, y::BigInt) = y*x

*(x::Integer, y::fmpz_mat) = BigInt(x)*y

*(x::fmpz_mat, y::Integer) = BigInt(y)*x

function +(x::fmpz_mat, y::Integer)
   z = parent(x)(x)
   for i = 1:min(rows(x), cols(x))
      z[i, i] += y
   end
   return z
end

+(x::Integer, y::fmpz_mat) = y + x

function -(x::fmpz_mat, y::Integer)
   z = parent(x)(x)
   for i = 1:min(rows(x), cols(x))
      z[i, i] -= y
   end
   return z
end

function -(x::Integer, y::fmpz_mat)
   z = parent(y)(y)
   for i = 1:min(rows(y), cols(y))
      z[i, i] = x - z[i, i]
   end
   return z
end

###########################################################################################
#
#   Scaling
#
###########################################################################################

function <<(x::fmpz_mat, y::Int)
   y < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpz_mat_scalar_mul_2exp, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Int), 
               &z, &x, y)
   return z
end

function >>(x::fmpz_mat, y::Int)
   y < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpz_mat_scalar_tdiv_q_2exp, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Int), 
               &z, &x, y)
   return z
end

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^(x::fmpz_mat, y::Int)
   y < 0 && throw(DomainError())
   rows(x) != cols(x) && error("Incompatible matrix dimensions")
   z = parent(x)()
   ccall((:fmpz_mat_pow, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Int), 
               &z, &x, y)
   return z
end

###########################################################################################
#
#   Comparisons
#
###########################################################################################

==(x::fmpz_mat, y::fmpz_mat) = ccall((:fmpz_mat_equal, :libflint), Bool, 
                                       (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &x, &y)

###########################################################################################
#
#   Ad hoc comparisons
#
###########################################################################################

function ==(x::fmpz_mat, y::Integer) 
   for i = 1:min(rows(x), cols(x))
      if x[i, i] != y
         return false
      end
   end
   for i = 1:rows(x)
      for j = 1:cols(x)
         if i != j && x[i, j] != 0
            return false
         end
      end
   end
   return true
end

==(x::Integer, y::fmpz_mat) = y == x

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv(x::fmpz_mat)
   z = parent(x)()
   d = fmpz()
   ccall((:fmpz_mat_inv, :libflint), Void, 
         (Ptr{fmpz_mat}, Ptr{fmpz}, Ptr{fmpz_mat}), &z, &d, &x)
   d1 = BigInt(d)
   if d1 == 1
      return z
   end
   if d1 == -1
      return -z
   end
   error("Matrix not invertible")
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

divexact(x::fmpz_mat, y::fmpz_mat) = x*inv(y)

###########################################################################################
#
#   Ad hoc exact division
#
###########################################################################################

function divexact(x::fmpz_mat, y::Int)
   z = parent(x)()
   ccall((:fmpz_mat_scalar_divexact_si, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Int), &z, &x, y)
   return z
end

function divexact(x::fmpz_mat, y::BigInt)
   z = parent(x)()
   temp = fmpz_readonly(y)
   ccall((:fmpz_mat_scalar_divexact_fmpz, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{fmpz_readonly}), &z, &x, &temp)
   return z
end

divexact(x::fmpz_mat, y::Integer) = divexact(x, BigInt(y))

###########################################################################################
#
#   Modular reduction
#
###########################################################################################

function reduce_mod(x::fmpz_mat, y::BigInt)
   z = parent(x)()
   temp = fmpz_readonly(y)
   ccall((:fmpz_mat_scalar_mod_fmpz, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{fmpz_readonly}), &z, &x, &temp)
   return z
end

reduce_mod(x::fmpz_mat, y::Integer) = reduce_mod(x, BigInt(y))

###########################################################################################
#
#   Characteristic polynomial
#
###########################################################################################

function charpoly(R::FmpzPolyRing, x::fmpz_mat)
   rows(x) != cols(x) && error("Non-square")
   z = R()
   ccall((:fmpz_mat_charpoly, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_mat}), &z, &x)
   return z
end

###########################################################################################
#
#   Determinant
#
###########################################################################################

function determinant(x::fmpz_mat)
   rows(x) != cols(x) && error("Non-square matrix")
   z = fmpz()
   ccall((:fmpz_mat_det, :libflint), Void, 
                (Ptr{fmpz}, Ptr{fmpz_mat}), &z, &x)
   return BigInt(z)
end

function determinant_divisor(x::fmpz_mat)
   rows(x) != cols(x) && error("Non-square matrix")
   z = fmpz()
   ccall((:fmpz_mat_det_divisor, :libflint), Void, 
                (Ptr{fmpz}, Ptr{fmpz_mat}), &z, &x)
   return BigInt(z)
end

function determinant_given_divisor(x::fmpz_mat, d::BigInt, proved=true)
   rows(x) != cols(x) && error("Non-square")
   z = fmpz()
   d1 = fmpz_readonly(d)
   ccall((:fmpz_mat_det_modular_given_divisor, :libflint), Void, 
                (Ptr{fmpz}, Ptr{fmpz_mat}, Ptr{fmpz_readonly}, Cint), &z, &x, &d1, proved)
   return BigInt(z)
end

function determinant_given_divisor(x::fmpz_mat, d::Integer, proved=true)
   return determinant_given_divisor(x, BigInt(d), proved)
end

###########################################################################################
#
#   Gram matrix
#
###########################################################################################

function gram(x::fmpz_mat)
   if rows(x) == cols(x)
      parz = parent(x)
   else
      parz = FmpzMatSpace(rows(x), rows(x))
   end
   z = parz()   
   ccall((:fmpz_mat_gram, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &x)
   return z
end

###########################################################################################
#
#   Hadamard matrix
#
###########################################################################################

function hadamard(R::FmpzMatSpace)
   R.rows != R.cols && error("Unable to create Hadamard matrix")
   z = R()
   success = ccall((:fmpz_mat_hadamard, :libflint), Bool, 
                   (Ptr{fmpz_mat},), &z)
   !success && error("Unable to create Hadamard matrix")
   return z
end

function is_hadamard(x::fmpz_mat)
   return ccall((:fmpz_mat_is_hadamard, :libflint), Bool, 
                   (Ptr{fmpz_mat},), &x)
end

###########################################################################################
#
#   Hermite normal form
#
###########################################################################################

function hnf(x::fmpz_mat)
   z = parent(x)()
   ccall((:fmpz_mat_hnf, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &x)
   return z
end

function hnf_with_transform(x::fmpz_mat)
   z = parent(x)()
   if rows(x) == cols(x)
      parz = parent(x)
   else
      parz = FmpzMatSpace(rows(x), rows(x))
   end
   u = parz()
   ccall((:fmpz_mat_hnf_transform, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &u, &x)
   return z, u
end

function hnf_modular(x::fmpz_mat, d::BigInt)
   z = parent(x)()
   d1 = fmpz_readonly(d)
   ccall((:fmpz_mat_hnf_modular, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{fmpz_readonly}), &z, &x, &d1)
   return z
end

function is_hnf(x::fmpz_mat)
   return ccall((:fmpz_mat_is_in_hnf, :libflint), Bool, 
                   (Ptr{fmpz_mat},), &x)
end

###########################################################################################
#
#   LLL
#
###########################################################################################

type lll_ctx
   delta::Float64
   eta::Float64
   rep_type::Int
   gram_type::Int

   function lll_ctx(delta::Float64, eta::Float64, rep=:zbasis, gram=:approx) 
      rt = rep == :zbasis ? 1 : 0
      gt = gram == :approx ? 0 : 1
      return new(delta, eta, rt, gt)
   end
end


function lll_with_transform(x::fmpz_mat, ctx=lll_ctx(0.99, 0.51))
   z = parent(x)(x)
   if rows(x) == cols(x)
      parz = parent(x)
   else
      parz = FmpzMatSpace(rows(x), rows(x))
   end
   u = parz(1)
   ccall((:fmpz_lll, :libflint), Void, 
         (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{lll_ctx}), &z, &u, &ctx)
   return z, u
end

function lll(x::fmpz_mat, ctx=lll_ctx(0.99, 0.51))
   z = parent(x)(x)
   if rows(x) == cols(x)
      parz = parent(x)
   else
      parz = FmpzMatSpace(rows(x), rows(x))
   end
   u = parz(1)
   ccall((:fmpz_lll, :libflint), Void, 
         (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{lll_ctx}), &z, &u, &ctx)
   return z
end

function lll_with_removal_transform(x::fmpz_mat, b::BigInt, ctx=lll_ctx(0.99, 0.51))
   z = parent(x)(x)
   if rows(x) == cols(x)
      parz = parent(x)
   else
      parz = FmpzMatSpace(rows(x), rows(x))
   end
   u = parz(1)
   b1 = fmpz_readonly(b)
   d = Int(ccall((:fmpz_lll_with_removal, :libflint), Cint, 
      (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{fmpz_readonly}, Ptr{lll_ctx}), &z, &u, &b1, &ctx))
   return d, z, u
end

function lll_with_removal(x::fmpz_mat, b::BigInt, ctx=lll_ctx(0.99, 0.51))
   z = parent(x)(x)
   if rows(x) == cols(x)
      parz = parent(x)
   else
      parz = FmpzMatSpace(rows(x), rows(x))
   end
   u = parz(1)
   b1 = fmpz_readonly(b)
   d = Int(ccall((:fmpz_lll_with_removal, :libflint), Cint, 
      (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{fmpz_readonly}, Ptr{lll_ctx}), &z, &u, &b1, &ctx))
   return d, z
end

###########################################################################################
#
#   Nullspace
#
###########################################################################################

function nullspace(x::fmpz_mat)
   z = parent(x)()
   if rows(x) == cols(x)
      parz = parent(x)
   else
      parz = FmpzMatSpace(cols(x), cols(x))
   end
   u = parz()
   rank = ccall((:fmpz_mat_nullspace, :libflint), Cint, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &u, &x)
   return u, rank
end

###########################################################################################
#
#   Rank
#
###########################################################################################

function rank(x::fmpz_mat)
   return ccall((:fmpz_mat_rank, :libflint), Int, 
                (Ptr{fmpz_mat},), &x)
end

###########################################################################################
#
#   Reduced row echelon form
#
###########################################################################################

function rref(x::fmpz_mat)
   z = parent(x)()
   d = fmpz()
   ccall((:fmpz_mat_rref, :libflint), Void,
         (Ptr{fmpz_mat}, Ptr{fmpz}, Ptr{fmpz_mat}), &z, &d, &x)
   return z, BigInt(d)
end

###########################################################################################
#
#   Smith normal form
#
###########################################################################################

function snf(x::fmpz_mat)
   z = parent(x)()
   ccall((:fmpz_mat_snf, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &x)
   return z
end

function snf_diagonal(x::fmpz_mat)
   z = parent(x)()
   ccall((:fmpz_mat_snf_diagonal, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &x)
   return z
end

function is_snf(x::fmpz_mat)
   return ccall((:fmpz_mat_is_in_snf, :libflint), Bool, 
                   (Ptr{fmpz_mat},), &x)
end

###########################################################################################
#
#   Linear solving
#
###########################################################################################

function solve(a::fmpz_mat, b::fmpz_mat)
   rows(a) != cols(a) && error("Not a square matrix in solve")
   rows(b) != rows(a) || cols(b) != 1 && ("Not a column vector in solve")
   z = parent(b)()
   d = fmpz()
   nonsing = ccall((:fmpz_mat_solve, :libflint), Bool,
         (Ptr{fmpz_mat}, Ptr{fmpz}, Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &d, &a, &b)
   !nonsing && error("Singular matrix in solve")
   return z, BigInt(d)
end

function solve_dixon(a::fmpz_mat, b::fmpz_mat)
   rows(a) != cols(a) && error("Not a square matrix in solve")
   rows(b) != rows(a) || cols(b) != 1 && ("Not a column vector in solve")
   z = parent(b)()
   d = fmpz()
   nonsing = ccall((:fmpz_mat_solve_dixon, :libflint), Bool,
         (Ptr{fmpz_mat}, Ptr{fmpz}, Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &d, &a, &b)
   !nonsing && error("Singular matrix in solve")
   return z, BigInt(d)
end

###########################################################################################
#
#   Trace
#
###########################################################################################

function trace(x::fmpz_mat)
   d = fmpz()
   return ccall((:fmpz_mat_trace, :libflint), Int, 
                (Ptr{fmpz}, Ptr{fmpz_mat}), &d, &x)
   return BigInt(d)
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function mul!(z::fmpz_mat, x::fmpz_mat, y::fmpz_mat)
   ccall((:fmpz_mat_mul, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &x, &y)
end

function addeq!(z::fmpz_mat, x::fmpz_mat)
   ccall((:fmpz_mat_add, :libflint), Void, 
                (Ptr{fmpz_mat}, Ptr{fmpz_mat}, Ptr{fmpz_mat}), &z, &z, &x)
end

###########################################################################################
#
#   Parent object call overloads
#
###########################################################################################

function Base.call(a::FmpzMatSpace)
   z = fmpz_mat(a.rows, a.cols)
   z.parent = a
   return z
end

function Base.call(a::FmpzMatSpace, arr::Array{BigInt, 2})
   z = fmpz_mat(a.rows, a.cols, arr)
   z.parent = a
   return z
end

function Base.call{T <: Integer}(a::FmpzMatSpace, arr::Array{T, 2})
   z = fmpz_mat(a.rows, a.cols, arr)
   z.parent = a
   return z
end

function Base.call(a::FmpzMatSpace, d::BigInt)
   z = fmpz_mat(a.rows, a.cols, d)
   z.parent = a
   return z
end

function Base.call(a::FmpzMatSpace, d::Integer)
   z = fmpz_mat(a.rows, a.cols, BigInt(d))
   z.parent = a
   return z
end

function Base.call(a::FmpzMatSpace, d::fmpz_mat)
   (a.rows != rows(d) || a.cols != cols(d)) && error("Incompatible matrix dimensions")
   z = fmpz_mat(d)
   z.parent = a
   return z
end

###########################################################################################
#
#   Promotions
#
###########################################################################################

Base.promote_rule{T <: Integer}(::Type{fmpz_mat}, ::Type{T}) = fmpz_mat

###########################################################################################
#
#   MatrixSpace constructor
#
###########################################################################################

function MatrixSpace(R::IntegerRing, r::Int, c::Int)
   return FmpzMatSpace(r, c)
end