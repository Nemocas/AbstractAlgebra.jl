###############################################################################
#
#   Matrix.jl : Generic mxn matrices over rings
#
###############################################################################

export Mat, MatrixSpace

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

elem_type{T <: RingElem}(::MatrixSpace{T}) = Mat{T}

base_ring(a::MatrixSpace) = a.base_ring

base_ring(a::MatElem) = base_ring(parent(a))

parent(a::MatElem) = a.parent

function check_parent(a::MatElem, b::MatElem)
   parent(a) != parent(b) && 
                error("Incompatible matrix spaces in matrix operation")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################    

function hash(a::MatElem)
   h = 0x3e4ea81eb31d94f4
   for i in 1:rows(a)
      for j in 1:cols(a)
         h $= hash(a[i, j])
         h = (h << 1) | (h >> (sizeof(Int)*8 - 1))
      end
   end
   return h
end

rows(a::MatElem) = parent(a).rows

cols(a::MatElem) = parent(a).cols

function getindex{T <: RingElem}(a::Mat{T}, r::Int, c::Int)
   return a.entries[r, c]
end
 
function setindex!{T <: RingElem}(a::Mat{T}, d::T, r::Int, c::Int)
   a.entries[r, c] = d
end

zero(a::MatrixSpace) = a()

one(a::MatrixSpace) = a(1)

function iszero(a::MatElem)
   for i = 1:rows(a)
      for j = 1:cols(a)
         if !iszero(a[i, j])
            return false
         end
      end
  end
  return true
end

function isone(a::MatElem)
   for i = 1:rows(a)
      for j = 1:cols(a)
         if i == j
            if !isone(a[i, j])
               return false
            end
         else
            if !iszero(a[i, j])
               return false
            end
         end
      end
  end
  return true
end

function deepcopy{T <: RingElem}(d::Mat{T})
   entries = Array(T, rows(d), cols(d))
   for i = 1:rows(d)
      for j = 1:cols(d)
         entries[i, j] = deepcopy(d[i, j])
      end
   end
   return parent(d)(entries)
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::MatElem) = canonical_unit(a[1, 1])

###############################################################################
#
#   AbstractString{} I/O
#
###############################################################################

function show(io::IO, a::MatrixSpace)
   print(io, "Matrix Space of ")
   print(io, a.rows, " rows and ", a.cols, " columns over ")
   print(io, base_ring(a))
end

function show(io::IO, a::MatElem)
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
         println(io, "")
      end
   end
end

show_minus_one{T <: RingElem}(::Type{Mat{T}}) = false

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(x::MatElem)
   par = parent(x)
   return par(-x.entries)
end

function transpose(x::MatElem)
   par = parent(x)
   return par(x.entries')
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +{T <: RingElem}(x::Mat{T}, y::Mat{T})
   check_parent(x, y)
   parz = parent(x)
   return parz(x.entries + y.entries)
end

function -{T <: RingElem}(x::Mat{T}, y::Mat{T})
   check_parent(x, y)
   parz = parent(x)
   return parz(x.entries - y.entries)
end

function *{T <: RingElem}(x::Mat{T}, y::Mat{T})
   cols(x) != rows(y) && error("Incompatible matrix dimensions")
   if rows(x) == cols(y) && rows(x) == cols(x)
      parz = parent(x)
   else
      parz = MatrixSpace(base_ring(x), rows(x), cols(y))
   end
   A = Array(T, rows(x), cols(y))
   C = base_ring(x)()
   for i = 1:rows(x)
      for j = 1:cols(y)
         A[i, j] = base_ring(x)()
         for k = 1:cols(x)
            mul!(C, x[i, k], y[k, j])
            addeq!(A[i, j], C)
         end
      end
   end
   return parz(A)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *{T <: RingElem}(x::Integer, y::Mat{T})
   z = similar(y.entries)
   parz = parent(y)
   for i = 1:rows(y)
      for j = 1:cols(y)
         z[i, j] = x*y[i, j]
      end
   end
   return parz(z)
end

function *{T <: RingElem}(x::fmpz, y::Mat{T})
   z = similar(y.entries)
   parz = parent(y)
   for i = 1:rows(y)
      for j = 1:cols(y)
         z[i, j] = x*y[i, j]
      end
   end
   return parz(z)
end

function *{T <: RingElem}(x::T, y::Mat{T})
   z = similar(y.entries)
   parz = parent(y)
   for i = 1:rows(y)
      for j = 1:cols(y)
         z[i, j] = x*y[i, j]
      end
   end
   return parz(z)
end

*{T <: RingElem}(x::Mat{T}, y::Integer) = y*x

*{T <: RingElem}(x::Mat{T}, y::fmpz) = y*x

*{T <: RingElem}(x::Mat{T}, y::T) = y*x

function +{T <: RingElem}(x::Integer, y::Mat{T})
   z = similar(y.entries)
   parz = parent(y)
   R = base_ring(y)
   for i = 1:rows(y)
      for j = 1:cols(y)
         if i != j
            z[i, j] = deepcopy(y[i, j])
         else
            z[i, j] = y[i, j] + R(x)
         end
      end
   end
   return parz(z)
end

+{T <: RingElem}(x::Mat{T}, y::Integer) = y + x

function +{T <: RingElem}(x::fmpz, y::Mat{T})
   z = similar(y.entries)
   parz = parent(y)
   R = base_ring(y)
   for i = 1:rows(y)
      for j = 1:cols(y)
         if i != j
            z[i, j] = deepcopy(y[i, j])
         else
            z[i, j] = y[i, j] + R(x)
         end
      end
   end
   return parz(z)
end

+{T <: RingElem}(x::Mat{T}, y::fmpz) = y + x

function +{T <: RingElem}(x::T, y::Mat{T})
   z = similar(y.entries)
   parz = parent(y)
   for i = 1:rows(y)
      for j = 1:cols(y)
         if i != j
            z[i, j] = deepcopy(y[i, j])
         else
            z[i, j] = y[i, j] + x
         end
      end
   end
   return parz(z)
end

+{T <: RingElem}(x::Mat{T}, y::T) = y + x

function -{T <: RingElem}(x::Integer, y::Mat{T})
   z = similar(y.entries)
   parz = parent(y)
   R = base_ring(y)
   for i = 1:rows(y)
      for j = 1:cols(y)
         if i != j
            z[i, j] = -y[i, j]
         else
            z[i, j] = R(x) - y[i, j] 
         end
      end
   end
   return parz(z)
end

function -{T <: RingElem}(x::Mat{T}, y::Integer) 
   z = similar(x.entries)
   parz = parent(x)
   R = base_ring(x)
   for i = 1:rows(x)
      for j = 1:cols(x)
         if i != j
            z[i, j] = deepcopy(x[i, j])
         else
            z[i, j] = x[i, j] - R(y)
         end
      end
   end
   return parz(z)
end

function -{T <: RingElem}(x::fmpz, y::Mat{T})
   z = similar(y.entries)
   parz = parent(y)
   R = base_ring(y)
   for i = 1:rows(y)
      for j = 1:cols(y)
         if i != j
            z[i, j] = -y[i, j]
         else
            z[i, j] = R(x) - y[i, j] 
         end
      end
   end
   return parz(z)
end

function -{T <: RingElem}(x::Mat{T}, y::fmpz) 
   z = similar(x.entries)
   parz = parent(x)
   R = base_ring(x)
   for i = 1:rows(x)
      for j = 1:cols(x)
         if i != j
            z[i, j] = deepcopy(x[i, j])
         else
            z[i, j] = x[i, j] - R(y)
         end
      end
   end
   return parz(z)
end

function -{T <: RingElem}(x::T, y::Mat{T})
   z = similar(y.entries)
   parz = parent(y)
   R = base_ring(y)
   for i = 1:rows(y)
      for j = 1:cols(y)
         if i != j
            z[i, j] = -y[i, j]
         else
            z[i, j] = x - y[i, j] 
         end
      end
   end
   return parz(z)
end

function -{T <: RingElem}(x::Mat{T}, y::T) 
   z = similar(x.entries)
   parz = parent(x)
   R = base_ring(x)
   for i = 1:rows(x)
      for j = 1:cols(x)
         if i != j
            z[i, j] = deepcopy(x[i, j])
         else
            z[i, j] = x[i, j] - y
         end
      end
   end
   return parz(z)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^{T <: RingElem}(a::Mat{T}, b::Int)
   b < 0 && throw(DomainError())
   rows(a) != cols(a) && error("Incompatible matrix dimensions in power")
   # special case powers of x for constructing polynomials efficiently
   if b == 0
      return one(parent(a))
   elseif b == 1
      return deepcopy(a)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit != 0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function =={T <: RingElem}(x::Mat{T}, y::Mat{T})
   check_parent(x, y)
   for i = 1:rows(x)
      for j = 1:cols(x)
         if x[i, j] != y[i, j]
            return false
         end
      end
   end
   return true
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function =={T <: RingElem}(x::Mat{T}, y::Integer) 
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

=={T <: RingElem}(x::Integer, y::Mat{T}) = y == x

function =={T <: RingElem}(x::Mat{T}, y::fmpz) 
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

=={T <: RingElem}(x::fmpz, y::Mat{T}) = y == x

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact{T <: RingElem}(x::Mat{T}, y::Integer)
   z = similar(x.entries)
   parz = parent(x)
   for i = 1:rows(x)
      for j = 1:cols(x)
         z[i, j] = divexact(x[i, j], y)
      end
   end
   return parz(z)
end

function divexact{T <: RingElem}(x::Mat{T}, y::fmpz)
   z = similar(x.entries)
   parz = parent(x)
   for i = 1:rows(x)
      for j = 1:cols(x)
         z[i, j] = divexact(x[i, j], y)
      end
   end
   return parz(z)
end

function divexact{T <: RingElem}(x::Mat{T}, y::T)
   z = similar(x.entries)
   parz = parent(x)
   for i = 1:rows(x)
      for j = 1:cols(x)
         z[i, j] = divexact(x[i, j], y)
      end
   end
   return parz(z)
end

###############################################################################
#
#   Gram matrix
#
###############################################################################

function gram{T <: RingElem}(x::Mat{T})
   if rows(x) == cols(x)
      parz = parent(x)
   else
      parz = MatrixSpace(base_ring(x), rows(x), rows(x))
   end
   z = parz()   
   for i = 1:rows(x)
      for j = 1:rows(x)
         z[i, j] = zero(base_ring(x))
         for k = 1:cols(x)
            z[i, j] += x[i, k]*x[j, k]
         end
      end
   end
   return z
end

###############################################################################
#
#   Trace
#
###############################################################################

function trace{T <: RingElem}(x::Mat{T})
   rows(x) != cols(x) && error("Not a square matrix in trace")
   d = zero(base_ring(x))
   for i = 1:rows(x)
      d += x[i, i]
   end
   return d
end

###############################################################################
#
#   Content
#
###############################################################################

function content{T <: RingElem}(x::Mat{T})
  d = zero(base_ring(x))
  for i = 1:rows(x)
     for j = 1:cols(x)
        d = gcd(d, x[i, j])
        if isone(d)
           return d
        end
     end
  end
  return d
end

###############################################################################
#
#   Permutation
#
###############################################################################

function *{T <: RingElem}(P::perm, x::Mat{T})
   z = parent(x)()
   m = rows(x)
   n = cols(x)
   for i = 1:m
      for j = 1:n
         z[P[i], j] = x[i, j]
      end
   end
   return z
end

###############################################################################
#
#   LU factorisation
#
###############################################################################

function lufact!{T <: FieldElem}(P::perm, A::Mat{T})
   m = rows(A)
   n = cols(A)
   rank = 0
   r = 1
   c = 1
   R = base_ring(A)
   t = R()
   while r <= m && c <= n
      if A[r, c] == 0
         i = r + 1
         while i <= m
            if A[i, c] != 0
               for j = 1:n
                  A.entries[i, j], A.entries[r, j] = A.entries[r, j], A.entries[i, j]
               end
               P[r], P[i] = P[i], P[r]
               break
            end
            i += 1
         end
         if i > m
            c += 1
            continue
         end
      end
      rank += 1
      d = -inv(A[r, c])
      for i = r + 1:m
         q = A[i, c]*d
         for j = c + 1:n
            mul!(t, A.entries[r, j], q)
            addeq!(A.entries[i, j], t)
         end
         A[i, c] = R()
         A[i, rank] = -q
      end
      r += 1
      c += 1
   end
   return rank
end

function lufact{T <: FieldElem}(P::perm, A::Mat{T})
   parent(P).n != rows(A) && error("Permutation does not match matrix")
   S = parent(A)
   R = base_ring(A)
   U = deepcopy(A)
   m = rows(A)
   n = cols(A)
   L = S()
   rank = lufact!(P, U)
   for i = 1:m
      for j = 1:n
         if i > j
            L[i, j] = U[i, j]
            U[i, j] = R()
         elseif i == j
            L[i, j] = R(1)
         else
            L[i, j] = R()
         end
      end
   end
   return rank, L, U
end

###############################################################################
#
#   Determinant
#
###############################################################################

function determinant_clow{T <: RingElem}(M::Mat{T})
   rows(M) != cols(M) && error("Dimensions don't match in determinant")
   R = base_ring(M)
   n = rows(M)
   A = Array(T, n, n)
   B = Array(T, n, n)
   C = R()
   for i = 1:n
      for j = 1:n
         A[i, j] = i == j ? R(1) : R(0)
         B[i, j] = R()
      end
   end
   for k = 1:n-1
      for i = 1:n
         for j = 1:i
            if !iszero(A[i, j])
               for m = j + 1:n
                  mul!(C, A[i, j], M[i, m])
                  addeq!(B[m, j], C)
               end
               for m = j + 1:n
                  mul!(C, A[i, j], M[i, j])
                  addeq!(B[m, m], -C)
               end
            end
         end
      end
      Temp = A
      A = B
      B = Temp
      if k != n - 1
         for i = 1:n
            for j = 1:i
               B[i, j] = R()
            end
         end
      end
   end
   D = R()
   for i = 1:n
      for j = 1:i
         if !iszero(A[i, j])
            D -= A[i, j]*M[i, j]
         end
      end
   end
   return isodd(n) ? -D : D
end

function charpoly{T <: RingElem}(V::PolynomialRing{T}, Y::Mat{T})
   rows(Y) != cols(Y) && error("Dimensions don't match in determinant")
   R = base_ring(Y)
   base_ring(V) != base_ring(Y) && error("Cannot coerce into polynomial ring")
   n = rows(Y)
   if n == 0
      return V()
   end
   F = Array(elem_type(R), n)
   A = Array(elem_type(R), n)
   M = Array(elem_type(R), n - 1, n)
   F[1] = -Y[1, 1]
   for i = 2:n
      F[i] = R()
      for j = 1:i
         M[1, j] = Y[j, i]
      end
      A[1] = Y[i, i]
      p = R()
      for j = 2:i - 1
         for k = 1:i
            s = R()
            for l = 1:i
               mul!(p, Y[k, l], M[j - 1, l])
               addeq!(s, p)
            end
            M[j, k] = s
         end
         A[j] = M[j, i]
      end 
      s = R()
      for j = 1:i
         mul!(p, Y[i, j], M[i - 1, j])
         addeq!(s, p)
      end
      A[i] = s
      for j = 1:i
         s = -F[j]
         for k = 1:j - 1
            mul!(p, A[k], F[j - k])
            addeq!(s, p)
         end
         F[j] = -s - A[j]
     end
   end
   z = gen(V)
   f = z^n
   for i = 1:n
      setcoeff!(f, n - i, F[i])
   end
   return f
end

function determinant{T <: RingElem}(M::Mat{T})
   R = base_ring(M)
   S, z = PolynomialRing(R, "z")
   n = rows(M)
   p = charpoly(S, M)
   d = coeff(p, 0)
   return isodd(n) ? -d : d
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

Base.promote_rule{T <: RingElem, V <: Integer}(::Type{Mat{T}}, ::Type{V}) = Mat{T}

Base.promote_rule{T <: RingElem}(::Type{Mat{T}}, ::Type{T}) = Mat{T}

Base.promote_rule{T <: RingElem}(::Type{Mat{T}}, ::Type{fmpz}) = Mat{T}

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call{T <: RingElem}(a::MatrixSpace{T}, b::RingElem)
   return a(base_ring(a)(b))
end

function Base.call{T <: RingElem}(a::MatrixSpace{T})
   entries = Array(T, a.rows, a.cols)
   for i = 1:a.rows
      for j = 1:a.cols
         entries[i, j] = zero(base_ring(a))
      end
   end
   z = Mat{T}(entries)
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::MatrixSpace{T}, b::Integer)
   entries = Array(T, a.rows, a.cols)
   for i = 1:a.rows
      for j = 1:a.cols
         if i != j
            entries[i, j] = zero(base_ring(a))
         else
            entries[i, j] = base_ring(a)(b)
         end
      end
   end
   z = Mat{T}(entries)
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::MatrixSpace{T}, b::T)
   parent(b) != base_ring(a) && error("Unable to coerce to matrix")
   entries = Array(T, a.rows, a.cols)
   for i = 1:a.rows
      for j = 1:a.cols
         if i != j
            entries[i, j] = zero(base_ring(a))
         else
            entries[i, j] = deepcopy(b)
         end
      end
   end
   z = Mat{T}(entries)
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::MatrixSpace{T}, b::Mat{T})
   parent(b) != a && error("Unable to coerce matrix")
   return b
end

function Base.call{T <: RingElem}(a::MatrixSpace{T}, b::Array{T, 2})
   if length(b) > 0
      parent(b[1, 1]) != base_ring(a) && error("Unable to coerce to matrix")
   end
   z = Mat{T}(b)
   z.parent = a
   return z
end

###############################################################################
#
#   MatrixSpace constructor
#
###############################################################################

function MatrixSpace(R::Ring, r::Int, c::Int)
   T = elem_type(R)
   return MatrixSpace{T}(R, r, c)
end

