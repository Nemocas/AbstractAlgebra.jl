################################################################################
#
#  nmod_mat.jl: flint nmod_mat types in julia
#
################################################################################

export nmod_mat, NmodMatSpace

export @check_parent

export getindex, setindex!, deepcopy, rows, cols, parent, base_ring, zero, one,
       issquare, show, ==, transpose, transpose!, -, +, *, ^, rref, rref!,
       trace, determinant, rank, inv, solve, lufact, sub, window, hcat, vcat,
       Array, lift, lift!, MatrixSpace

################################################################################
#
#  Macro for parent checking
#
################################################################################

macro check_parent(x,y)
  return :(parent($x) == parent($y) ? nothing :
                      error("Arguments have different parents"))
end

################################################################################
#
#  Types and memory management
#
################################################################################

NmodMatID = ObjectIdDict()

type NmodMatSpace <: Ring
  base_ring::ResidueRing{fmpz}
  _n::UInt
  rows::Int
  cols::Int

  function NmodMatSpace(R::ResidueRing{fmpz}, r::Int, c::Int)
    (r <= 0 || c <= 0) && error("Dimensions must be positive")
    ZZ(typemax(UInt)) < abs(R.modulus) &&
      error("Modulus of ResidueRing must less then ", ZZ(typemax(UInt)))
    try
      return NmodMatID[R, r, c]
    catch
      NmodMatID[R, r, c] = new(R, UInt(BigInt(abs(R.modulus))), r, c)
    end
  end
end

type nmod_mat <: MatElem
  entries::Ptr{Void}
  r::Int                  # Clong
  c::Int                  # Clong
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
        ccall((:_nmod_mat_set_entry, :libflint), Void,
                (Ptr{nmod_mat}, Int, Int, UInt),
                &z, i-1, j-1, arr[i,j] % n)
      ## I am not sure if argument must be already reduced
      end
    end
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{fmpz, 2})
    (r <= 0 || c <= 0) && error("Dimensions must be positive")
    (size(arr) != (r,c)) && error("Array of wrong dimension")
    z = new()
    t = ZZ()
    tt = UInt(0)
    ccall((:nmod_mat_init, :libflint), Void,
            (Ptr{nmod_mat}, Int, Int, UInt), &z, r, c, n)
    finalizer(z, _nmod_mat_clear_fn)
    for i = 1:r
      for j = 1:c
        ccall((:fmpz_mod_ui, :libflint), UInt,
                (Ptr{fmpz}, Ptr{fmpz}, UInt), &t, &arr[i,j], n)
        tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
        ccall((:_nmod_mat_set_entry, :libflint), Void,
                (Ptr{nmod_mat}, Int, Int, UInt),
                &z, i-1, j-1, tt)
      end
    end
    return z
  end

  function nmod_mat{T <: Integer}(r::Int, c::Int, n::UInt, arr::Array{T, 2})
    arr = map(ZZ, arr)
    return nmod_mat(r, c, n, arr)
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{Residue{fmpz}, 2})
    (r <= 0 || c <= 0) && error("Dimensions must be positive")
    (size(arr) != (r,c)) && error("Array of wrong dimension")
    z = new()
    t = ZZ()
    tt = UInt(0)
    ccall((:nmod_mat_init, :libflint), Void,
            (Ptr{nmod_mat}, Int, Int, UInt), &z, r, c, n)
    finalizer(z, _nmod_mat_clear_fn)
    for i = 1:r
      for j = 1:c
        ccall((:fmpz_mod_ui, :libflint), UInt,
                (Ptr{fmpz}, Ptr{fmpz}, UInt), &t, &arr[i,j].data, n)
        tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
        ccall((:_nmod_mat_set_entry, :libflint), Void,
                (Ptr{nmod_mat}, Int, Int, UInt),
                &z, i-1, j-1, tt)
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
    (n > ZZ(typemax(UInt))) &&
          error("Exponent must be smaller then ", ZZ(typemax(UInt)))
    return nmod_mat(UInt(n), b) 
  end
end

function _nmod_mat_clear_fn(mat::nmod_mat)
  ccall((:nmod_mat_clear, :libflint), Void, (Ptr{nmod_mat}, ), &mat)
end

function getindex(a::nmod_mat, i::Int, j::Int)
  checkbounds(a.r, i)
  checkbounds(a.c, j)
  u = ccall((:_nmod_mat_get_entry, :libflint), UInt,
              (Ptr{nmod_mat}, Int, Int), &a, i - 1 , j - 1)
  return base_ring(a)(u)
end

function setindex!(a::nmod_mat, u::UInt, i::Int, j::Int)
  checkbounds(a.r, i)
  checkbounds(a.c, j)
  ccall((:_nmod_mat_set_entry, :libflint), Void,
          (Ptr{nmod_mat}, Int, Int, UInt), &a, i-1, j-1, u)
end

function setindex!(a::nmod_mat, u::fmpz, i::Int, j::Int)
  checkbounds(a.r, i)
  checkbounds(a.c, j)
  t = ZZ()
  ccall((:fmpz_mod_ui, :libflint), UInt,
          (Ptr{fmpz}, Ptr{fmpz}, UInt), &t, &u, a._n)
  tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
  setindex!(a,tt,i,j)
end

function setindex!(a::nmod_mat, u::Residue{fmpz}, i::Int, j::Int)
  checkbounds(a.r, i)
  checkbounds(a.c, j)
  (base_ring(a) != parent(u)) && error("Parent objects must coincide") 
  setindex!(a, u.data, i, j)
end

setindex!(a::nmod_mat, u::Integer, i::Int, j::Int) =
        setindex!(a, fmpz(u), i, j)


function deepcopy(a::nmod_mat)
  z = nmod_mat(rows(a), cols(a), a._n)
  if isdefined(a, :parent)
    z.parent = a.parent
  end
  ccall((:nmod_mat_set, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &a)
  return z
end

rows(a::nmod_mat) = a.r

cols(a::nmod_mat) = a.c

parent(a::nmod_mat) = a.parent

base_ring(a::NmodMatSpace) = a.base_ring

base_ring(a::nmod_mat) = a.parent.base_ring

zero(a::NmodMatSpace) = a()

function one(a::NmodMatSpace)
  (a.rows != a.cols) && error("Matrices must be quadratic")
  z = a()
  ccall((:nmod_mat_one, :libflint), Void, (Ptr{nmod_mat}, ), &z)
  return z
end

################################################################################
#
#  Predicates
#
################################################################################

issquare(a::nmod_mat) = (rows(a) == cols(a))

function iszero(a::nmod_mat)
  r = ccall((:nmod_mat_is_zero, :libflint), Cint, (Ptr{nmod_mat}, ), &a)
  return Bool(r)
end

################################################################################
#
#  String I/O
#
################################################################################

function show(io::IO, a::NmodMatSpace)
   print(io, "Matrix Space of ")
   print(io, a.rows, " rows and ", a.cols, " columns over ")
   print(io, a.base_ring)
end

function show(io::IO, a::nmod_mat)
   rows = a.r
   cols = a.c
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
         println(io, "")
      end
   end
end

################################################################################
#
#  Comparison
#
################################################################################

==(a::nmod_mat, b::nmod_mat) = (a.parent == b.parent) &&
        Bool(ccall((:nmod_mat_equal, :libflint), Cint,
                (Ptr{nmod_mat}, Ptr{nmod_mat}), &a, &b))

################################################################################
#
#  Transpose
#
################################################################################

# are we allowing non-square matrices?!

function transpose(a::nmod_mat)
  z = NmodMatSpace(base_ring(a), parent(a).cols, parent(a).rows)()
  ccall((:nmod_mat_transpose, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &a)
  return z
end

function transpose!(a::nmod_mat)
  !issquare(a) && error("Matrix must be a square matrix")
  ccall((:nmod_mat_transpose, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}), &a, &a)
end

################################################################################
#
#  Arithmetic
#
################################################################################

function -(x::nmod_mat)
  z = parent(x)()
  ccall((:nmod_mat_neg, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &x)
  return z
end

function +(x::nmod_mat, y::nmod_mat)
  @check_parent(x,y)
  z = parent(x)()
  ccall((:nmod_mat_add, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &x, &y)
  return z
end

function -(x::nmod_mat, y::nmod_mat)
  @check_parent(x,y)
  z = parent(x)()
  ccall((:nmod_mat_sub, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &x, &y)
  return z
end

function *(x::nmod_mat, y::nmod_mat)
  (base_ring(x) != base_ring(y)) && error("Base ring must be equal")
  (cols(x) != rows(y)) && error("Dimensions are wrong")
  z = MatrixSpace(base_ring(x), rows(x), cols(y))()
  ccall((:nmod_mat_mul, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &x, &y)
  return z
end

function *(x::nmod_mat, y::UInt)
  z = parent(x)()
  ccall((:nmod_mat_scalar_mul, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}, UInt), &z, &x, y)
  return z
end

*(x::UInt, y::nmod_mat) = y*x

function *(x::nmod_mat, y::fmpz)
  t = ZZ()
  ccall((:fmpz_mod_ui, :libflint), UInt,
          (Ptr{fmpz}, Ptr{fmpz}, UInt), &t, &y, parent(x)._n)
  tt = ccall((:fmpz_get_ui, :libflint), UInt, (Ptr{fmpz}, ), &t)
  return x*tt
end

*(x::fmpz, y::nmod_mat) = y*x

function *(x::nmod_mat, y::Integer)
  return x*fmpz(y)
end

*(x::Integer, y::nmod_mat) = y*x

function *(x::nmod_mat, y::Residue{fmpz})
  (base_ring(x) != parent(y)) && error("Parent objects must coincide")
  return x*y.data
end

*(x::Residue{fmpz}, y::nmod_mat) = y*x

function ^(x::nmod_mat, y::UInt)
  z = parent(x)()
  ccall((:nmod_mat_pow, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}, UInt), &z, &x, y)
  return z
end

function ^(x::nmod_mat, y::Int)
  ( y < 0 ) && error("Exponent must be positive")
  return x^UInt(y)
end

function ^(x::nmod_mat, y::fmpz)
  ( y < 0 ) && error("Exponent must be positive")
  ( y > ZZ(typemax(UInt))) &&
          error("Exponent must be smaller then ", ZZ(typemax(UInt)))
  return x^(UInt(y))
end

################################################################################
#
#  Row echelon form
#
################################################################################

function rref(a::nmod_mat)
  z = deepcopy(a)
  ccall((:nmod_mat_rref, :libflint), Void, (Ptr{nmod_mat}, ), &z)
  return z
end

function rref!(a::nmod_mat)
  ccall((:nmod_mat_rref, :libflint), Void, (Ptr{nmod_mat}, ), &a)
  return a
end

################################################################################
#
#  Trace and Determinant
#
################################################################################

function trace(a::nmod_mat)
  !issquare(a) && error("Matrix must be a square matrix")
  r = ccall((:nmod_mat_trace, :libflint), UInt, (Ptr{nmod_mat}, ), &a)
  return base_ring(a)(r)
end

function determinant(a::nmod_mat)
  !issquare(a) && error("Matrix must be a square matrix")
  r = ccall((:nmod_mat_det, :libflint), UInt, (Ptr{nmod_mat}, ), &a)
  return base_ring(a)(r)
end

################################################################################
#
#  Rank
#
################################################################################

function rank(a::nmod_mat)
  r = ccall((:nmod_mat_rank, :libflint), Int, (Ptr{nmod_mat}, ), &a)
  return r
end

################################################################################
#
#  Inverse
#
################################################################################

function inv(a::nmod_mat)
  !issquare(a) && error("Matrix must be a square matrix")
  z = parent(a)()
  r = ccall((:nmod_mat_inv, :libflint), Int,
          (Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &a)
  !Bool(r) && error("Matrix not invertible")
  return z
end

################################################################################
#
#  Linear solving
#
################################################################################

function solve(x::nmod_mat, y::nmod_mat)
  (base_ring(x) != base_ring(y)) && error("Matrices must have same base ring")
  !issquare(x)&& error("First argument not a square matrix in solve")
  (y.r != x.r) || y.c != 1 && ("Not a column vector in solve")
  z = parent(y)()
  r = ccall((:nmod_mat_solve, :libflint), Int,
          (Ptr{nmod_mat}, Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &x, &y)
  !Bool(r) && error("Singular matrix in solve")
  return z
end

################################################################################
#
#  LU decomposition
#
################################################################################

function lufact(x::nmod_mat)
  t = deepcopy(x)
  p = Array(Int, x.r)
  r = ccall((:nmod_mat_lu, :libflint), Cint,
          (Ptr{Int}, Ptr{nmod_mat}, Cint), p, &t, 0)
  r = Int(r)
  if issquare(x) && r == rows(x)
    l = deepcopy(t)
    for i in 1:cols(l)
      l[i,i] = 1
    end
    for i in 1:rows(l)
      for j in i+1:cols(l)
        l[i,j] = 0
      end
    end
    for i in 1:cols(t)
      for j in 1:i-1
        t[i,j] = 0
      end
    end
    u = t
  else
    l = window(t, 1, 1, rows(x), r)
    for i in 1:r 
      l[i,i] = 1
    end
    for i in 1:rows(l)
      for j in i+1:cols(l)
        l[i,j] = 0
      end
    end
    u = window(t, 1, 1, r, cols(x))
      for i in 1:rows(u)
        for j in 1:i-1
          u[i,j] = 0
        end
      end
    end
  return l,u,p
end

################################################################################
#
#  Windowing !!! Not documented, but useful, flint functions
#
################################################################################

function window(x::nmod_mat, r1::Int, c1::Int, r2::Int, c2::Int)
  checkbounds(x.r, r1)
  checkbounds(x.r, r2)
  checkbounds(x.c, c1)
  checkbounds(x.c, c2)
  (r1 > r2 || c1 > c2) && error("Invalid parameters")
  temp = MatrixSpace(parent(x).base_ring, r2 - r1 + 1, c2 - c1 + 1)()
  ccall((:nmod_mat_window_init, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}, Int, Int, Int, Int),
          &temp, &x, r1-1, c1-1, r2, c2)
  z = deepcopy(temp)
  ccall((:nmod_mat_window_clear, :libflint), Void, (Ptr{nmod_mat}, ), &temp)
  return z
end

function window(x::nmod_mat, r::UnitRange{Int}, c::UnitRange{Int})
  return window(x, r.start, c.start, r.stop, c.stop)
end

sub(x::nmod_mat, r1::Int, c1::Int, r2::Int, c2::Int) =
        window(x, r1, c1, r2, c2)

sub(x::nmod_mat, r::UnitRange{Int}, c::UnitRange{Int}) = window(x, r, c)
  
################################################################################
#
#  Concatenation !!! Not documented, but useful, flint functions
#
################################################################################

function hcat(x::nmod_mat, y::nmod_mat)
  (base_ring(x) != base_ring(y)) && error("Matrices must have same base ring")
  (x.r != y.r) && error("Matrices must have same number of rows")
  z = MatrixSpace(base_ring(x), x.r, x.c + y.c)()
  ccall((:nmod_mat_concat_horizontal, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &x, &y)
  return z
end

function vcat(x::nmod_mat, y::nmod_mat)
  (base_ring(x) != base_ring(y)) && error("Matrices must have same base ring")
  (x.c != y.c) && error("Matrices must have same number of columns")
  z = MatrixSpace(base_ring(x), x.r + y.r, x.c)()
  ccall((:nmod_mat_concat_vertical, :libflint), Void,
          (Ptr{nmod_mat}, Ptr{nmod_mat}, Ptr{nmod_mat}), &z, &x, &y)
  return z
end

################################################################################
#
#  Conversion
#
################################################################################

function Array(b::nmod_mat)
  a = Array(Residue{fmpz}, b.r, b.c)
  for i = 1:b.r
    for j = 1:b.c
      a[i,j] = base_ring(b)(ccall((:_nmod_mat_get_entry, :libflint), UInt,
              (Ptr{nmod_mat}, Int, Int), &b, i-1, j-1))
    end
  end
  return a
end

function lift(a::nmod_mat)
  z = MatrixSpace(ZZ, rows(a), cols(a))()
  ccall((:fmpz_mat_set_nmod_mat, :libflint), Void,
          (Ptr{fmpz_mat}, Ptr{nmod_mat}), &z, &a)
  return z 
end

function lift!(z::fmpz_mat, a::nmod_mat)
  ccall((:fmpz_mat_set_nmod_mat, :libflint), Void,
          (Ptr{fmpz_mat}, Ptr{nmod_mat}), &z, &a)
  return z 
end

################################################################################
#
#  Parent object overloading
#
################################################################################

function Base.call(a::NmodMatSpace)
  z = nmod_mat(a.rows, a.cols, a._n)
  z.parent = a
  return z
end

function Base.call(a::NmodMatSpace, arr::Array{BigInt, 2})
  z = nmod_mat(a.rows, a.cols, a._n, arr)
  z.parent = a
  return z
end

function Base.call(a::NmodMatSpace, arr::Array{fmpz, 2})
  z = nmod_mat(a.rows, a.cols, a._n, arr)
  z.parent = a
  return z
end

function Base.call(a::NmodMatSpace, arr::Array{Int, 2})
  z = nmod_mat(a.rows, a.cols, a._n, arr)
  z.parent = a
  return z
end

function Base.call(a::NmodMatSpace, arr::Array{Residue{fmpz}, 2})
  length(arr) == 0 && error("Array must be nonempty")
  (base_ring(a) != parent(arr[1])) && error("Elements must have same base ring")
  z = nmod_mat(a.rows, a.cols, a._n, arr)
  z.parent = a
  return z
end

function Base.call(a::NmodMatSpace, arr::Array{Int, 1})
  (length(arr) != a.cols * a.rows) &&
          error("Array must be of length ", a.cols * a.rows)

  arr = transpose(reshape(arr,a.cols,a.rows))
  return a(arr)
end

function Base.call(a::NmodMatSpace, arr::Array{BigInt, 1})
  (length(arr) != a.cols * a.rows) &&
          error("Array must be of length ", a.cols * a.rows)
  arr = transpose(reshape(arr,a.cols,a.rows))
  return a(arr)
end

function Base.call(a::NmodMatSpace, arr::Array{fmpz, 1})
  (length(arr) != a.cols * a.rows) &&
          error("Array must be of length ", a.cols * a.rows)
  arr = transpose(reshape(arr,a.cols,a.rows))
  return a(arr)
end

function Base.call(a::NmodMatSpace, arr::Array{Residue{fmpz}, 1})
  (length(arr) != a.cols * a.rows) &&
          error("Array must be of length ", a.cols * a.rows)
  arr = transpose(reshape(arr,a.cols,a.rows))
  return a(arr)
end

function Base.call(a::NmodMatSpace, b::fmpz_mat)
  (a.cols != b.c || a.rows != b.r) && error("Dimensions do not fit")
  z = nmod_mat(a._n, b)
  z.parent = a
  return z
end

################################################################################
#
#  Matrix space constructor
#
################################################################################

function MatrixSpace(R::ResidueRing{fmpz}, r::Int, c::Int)
  return try
    NmodMatSpace(R, r, c)
  catch
    error("Not yet implemented")
  end
end

################################################################################
#
#  Helper functions, should go into ZZ.jl
#
################################################################################

convert(::Type{UInt64}, x::fmpz) =
        ccall((:fmpz_get_ui, :libflint), UInt64, (Ptr{fmpz}, ), &x)
