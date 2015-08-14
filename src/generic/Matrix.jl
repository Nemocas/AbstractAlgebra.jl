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

