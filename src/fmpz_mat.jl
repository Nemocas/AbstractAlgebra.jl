###########################################################################################
#
#   fmpz_mat.jl : Flint matrices over ZZ
#
###########################################################################################

import Base: getindex, setindex!

export fmpz_mat, MatrixSpace, getindex, setindex!

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

   function FmpzMatSpace(r::Int, c::Int)
      return try
         FmpzMatID[r, c]
      catch
         FmpzMatID[r, c] = new(r, c)
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
end

function _fmpz_mat_clear_fn(a::fmpz_mat)
   ccall((:fmpz_mat_clear, :libflint), Void, (Ptr{fmpz_mat},), &a)
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function getindex(a::fmpz_mat, r::Int, c::Int)
   (r > a.parent.rows || c > a.parent.rows) && throw(BoundsError())
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

###########################################################################################
#
#   MatrixSpace constructor
#
###########################################################################################

function MatrixSpace(R::IntegerRing, r::Int, c::Int)
   return FmpzMatSpace(r, c)
end