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
#   String I/O
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
      parz = FmpzMatSpace(rows(x), cols(y))
   end
   return parz(x.entries*y.entries)
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
#   Promotion rules
#
###############################################################################

Base.promote_rule{T <: RingElem, V <: Integer}(::Type{Mat{T}}, ::Type{V}) = Mat{T}

Base.promote_rule{T <: RingElem}(::Type{Mat{T}}, ::Type{T}) = Mat{T}

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

