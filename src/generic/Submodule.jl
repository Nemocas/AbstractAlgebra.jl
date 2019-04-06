###############################################################################
#
#   Submodule.jl : Submodules of modules
#
###############################################################################

export Submodule, submodule_elem, ngens, supermodule

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{submodule_elem{T}}) where T <: RingElement = Submodule{T}

elem_type(::Type{Submodule{T}}) where T <: RingElement = submodule_elem{T}

parent(v::submodule_elem) = v.parent

base_ring(N::Submodule{T}) where T <: RingElement = N.base_ring

base_ring(v::submodule_elem{T}) where T <: RingElement = base_ring(v.parent)

ngens(N::Submodule{T}) where T <: RingElement = length(N.gens)

gens(N::Submodule{T}) where T <: RingElement = [gen(N, i) for i = 1:ngens(N)]

function gen(N::Submodule{T}, i::Int) where T <: RingElement
   R = base_ring(N)
   return N([(j == i ? one(R) : zero(R)) for j = 1:ngens(N)])
end

@doc Markdown.doc"""
    supermodule(M::Submodule{T}) where T <: RingElement
> Return the module that this module is a submodule of.
"""
supermodule(M::Submodule{T}) where T <: RingElement = M.m

function check_parent(v1::submodule_elem{T}, v2::submodule_elem{T}) where T <: RingElement
   parent(v1) != parent(v2) && error("Incompatible module elements")
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, N::Submodule{T}) where T <: RingElement
   println(io, "Submodule of:")
   print(IOContext(io, :compact => true), N.m)
   println(io, "")
   println(io, " with generators:")
   print(IOContext(io, :compact => true), N.gens)
end

function show(io::IO, N::Submodule{T}) where T <: FieldElement
   println(io, "Subspace of:")
   print(IOContext(io, :compact => true), N.m)
   println(io, "")
   println(io, " with generators:")
   print(IOContext(io, :compact => true), N.gens)
end

function show(io::IO, v::submodule_elem)
   print(io, "(")
   len = ngens(parent(v))
   for i = 1:len - 1
      print(IOContext(io, :compact => true), v.v[1, i])
      print(io, ", ")
   end
   if len > 0
      print(IOContext(io, :compact => true), v.v[1, len])
   end
   print(io, ")")
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(v::submodule_elem{T}) where T <: RingElement
   N = parent(v)
   return N(-v.v)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(v1::submodule_elem{T}, v2::submodule_elem{T}) where T <: RingElement
   check_parent(v1, v2)
   N = parent(v1)
   return N(v1.v + v2.v)
end

function -(v1::submodule_elem{T}, v2::submodule_elem{T}) where T <: RingElement
   check_parent(v1, v2)
   N = parent(v1)
   return N(v1.v - v2.v)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(v::submodule_elem{T}, c::T) where T <: RingElem
   N = parent(v)
   return N(v.v*c)
end

function *(v::submodule_elem{T}, c::U) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(v.v*c)
end

function *(c::T, v::submodule_elem{T}) where T <: RingElem
   N = parent(v)
   return N(c*v.v)
end

function *(c::U, v::submodule_elem{T}) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(c*v.v)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(m::submodule_elem{T}, n::submodule_elem{T}) where T <: RingElement
   check_parent(m, n)
   return m.v == n.v
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (N::Submodule{T})(v::Vector{T}) where T <: RingElement
   length(v) != ngens(N) && error("Length of vector does not match number of generators")
   mat = matrix(base_ring(N), 1, length(v), v)
   return submodule_elem{T}(N, mat)
end

function (N::Submodule{T})(v::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(v) != ngens(N) && error("Length of vector does not match number of generators")
   nrows(v) != 1 && ("Not a vector in submodule_elem constructor")
   return submodule_elem{T}(N, v)
end

###############################################################################
#
#   Submodule constructor
#
###############################################################################

reduced_form(mat::AbstractAlgebra.MatElem{T}) where T <: RingElement = hnf(mat)

function reduced_form(mat::AbstractAlgebra.MatElem{T}) where T <: FieldElement
  r, m = rref(mat)
  return mat
end

@doc Markdown.doc"""
    Submodule(m::AbstractAlgebra.Module{T}, gens::Vector{<:AbstractAlgebra.ModuleElem{T}}) where T <: RingElement
> Return the submodule of the module `m` generated by the given generators,
> given as elements of `m`.
"""
function Submodule(m::AbstractAlgebra.Module{T}, gens::Vector{<:AbstractAlgebra.ModuleElem{T}}) where T <: RingElement
   if length(gens) == 0
      M = Submodule{T}(m, gens)
      f = map_from_func(M, m, x -> nothing)
      return M, f
   end
   r = length(gens)
   s = ngens(m)
   mat = matrix(base_ring(m), r, s, [0 for i in 1:r*s])
   for i = 1:r
      parent(gens[i]) != m && error("Incompatible module elements")
      for j = 1:s
         mat[i, j] = gens[i].v[1, j]
      end
   end
   mat = reduced_form(mat)
   num = r
   while num > 0
      rowzero = true
      for j = 1:s
         if mat[num, j] != 0
            rowzero = false
         end
      end
      if !rowzero
         break
      end
      num -= 1
   end
   gens = [m([mat[i, j] for j = 1:s]) for i = 1:num]
   M = Submodule{T}(m, gens)
   f = map_from_func(M, m, x -> sum(x.v[1, i]*gens[i] for i in 1:ncols(x.v)))
   M.map = f
   return M, f
end

