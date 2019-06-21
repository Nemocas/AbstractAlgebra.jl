###############################################################################
#
#   Submodule.jl : Submodules of modules
#
###############################################################################

export Submodule, submodule_elem, ngens, gens, supermodule

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

ngens(N::Submodule{T}) where T <: RingElement = length(N.gen_cols)

gens(N::Submodule{T}) where T <: RingElement = [gen(N, i) for i = 1:ngens(N)]

function gen(N::Submodule{T}, i::Int) where T <: RingElement
   R = base_ring(N)
   return N([(j == i ? one(R) : zero(R)) for j = 1:ngens(N)])
end

# Generators as elements of supermodule. Used internally.
generators(N::Submodule{T}) where T <: RingElement = N.gens::Vector{elem_type(N.m)}

@doc Markdown.doc"""
    supermodule(M::Submodule{T}) where T <: RingElement
> Return the module that this module is a submodule of.
"""
supermodule(M::Submodule{T}) where T <: RingElement = M.m

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, N::Submodule{T}) where T <: RingElement
   print(io, "Submodule over ")
   print(IOContext(io, :compact => true), base_ring(N))
   show_gens_rels(io, N)
end

function show(io::IO, N::Submodule{T}) where T <: FieldElement
   print(io, "Subspace over ")
   print(IOContext(io, :compact => true), base_ring(N))
   show_gens_rels(io, N)
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
   base_ring(v) != parent(c) && error("Incompatible rings")
   N = parent(v)
   return N(v.v*c)
end

function *(v::submodule_elem{T}, c::U) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(v.v*c)
end

function *(c::T, v::submodule_elem{T}) where T <: RingElem
   base_ring(v) != parent(c) && error("Incompatible rings")
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
   mat = reduce_mod_rels(mat, rels(N), 1)
   return submodule_elem{T}(N, mat)
end

function (N::Submodule{T})(v::Vector{Any}) where T <: RingElement
   length(v) != 0 && error("Incompatible element")
   return N(T[])
end

function (N::Submodule{T})(v::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(v) != ngens(N) && error("Length of vector does not match number of generators")
   nrows(v) != 1 && ("Not a vector in submodule_elem constructor")
   v = reduce_mod_rels(v, rels(N), 1)
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
  return m
end

isreduced_form(mat::AbstractAlgebra.MatElem{T}) where T <: RingElement = ishnf(mat)

isreduced_form(mat::AbstractAlgebra.MatElem{T}) where T <: FieldElement = isrref(mat)

@doc Markdown.doc"""
    Submodule(m::AbstractAlgebra.FPModule{T}, gens::Vector{<:AbstractAlgebra.FPModuleElem{T}}) where T <: RingElement
> Return the submodule of the module `m` generated by the given generators,
> given as elements of `m`.
"""
function Submodule(m::AbstractAlgebra.FPModule{T}, gens::Vector{S}) where {T <: RingElement, S <: AbstractAlgebra.FPModuleElem{T}}
   R = base_ring(m)
   r = length(gens)
   while r > 0 && iszero_row(gens[r].v, 1) # check that not all gens are zero
      r -= 1
   end
   if r == 0
      gens = Vector{S}(undef, 0) # original may have generators that are zero
      M = Submodule{T}(m, gens, Vector{dense_matrix_type(T)}(undef, 0),
                       Vector{Int}(undef, 0), Vector{Int}(undef, 0))
      f = ModuleHomomorphism(M, m, matrix(R, 0, ngens(m), []))
      M.map = f
      return M, f
   end
   # Make generators rows of a matrix
   s = ngens(m)
   mat = zero_matrix(base_ring(m), r, s)
   for i = 1:r
      parent(gens[i]) !== m && error("Incompatible module elements")
      for j = 1:s
         mat[i, j] = gens[i].v[1, j]
      end
   end
   # Reduce matrix (hnf/rref), ensuring we remain reduced mod old relations
   old_rels = rels(m)
   while !isreduced_form(mat)
      # Reduce matrix (hnf/rref)
      mat = reduced_form(mat)
      # Remove zero rows
      num = r
      while num > 0 && iszero_row(mat, num)
         num -= 1
      end
      # Reduce modulo old relations
      for i = 1:num
         Mi = @view mat[i:i, :]
         g = reduce_mod_rels(Mi, old_rels, 1)
         for j = 1:ncols(Mi)
            Mi[1, j] = g[1, j]
         end
      end 
   end
   # Remove zero rows
   num = r
   while num > 0 && iszero_row(mat, num)
      num -= 1
   end
   # Rewrite matrix without zero rows and add old relations as rows
   # We flip the rows so the output of kernel is upper triangular with
   # respect to the original data, which saves time in reduced_form
   nr = num + length(old_rels)
   new_mat = zero_matrix(base_ring(m), nr, s)
   for i = 1:num
      for j = 1:s
         new_mat[nr - i + 1, j] = mat[i, j]
      end
   end
   for i = 1:length(old_rels)
      for j = 1:s
         new_mat[nr - i - num + 1, j] = old_rels[i][1, j]
      end
   end
   # Rewrite old relations in terms of generators of new submodule
   num_rels, K = left_kernel(new_mat)
   new_rels = zero_matrix(base_ring(m), num_rels, num)
   # we flip rows and columns so that input is in terms of original data and
   # in upper triangular form, to save time in reduced_form below
   for j = 1:num_rels
      for k = 1:num
         new_rels[num_rels - j + 1, k] = K[j, nr - k + 1]
      end
   end
   # Compute reduced form of new rels
   new_rels = reduced_form(new_rels)
   # remove rows and columns corresponding to unit pivots
   gen_cols, culled, pivots = cull_matrix(new_rels)
   # put all the culled relations into new relations
   T1 = typeof(new_rels)
   srels = Vector{T1}(undef, length(culled))
   for i = 1:length(culled)
      srels[i] = matrix(R, 1, length(gen_cols),
                    [new_rels[culled[i], gen_cols[j]]
                       for j in 1:length(gen_cols)])
   end
   # Make submodule whose generators are the nonzero rows of mat
   nonzero_gens = [m([mat[i, j] for j = 1:s]) for i = 1:num]
   M = Submodule{T}(m, nonzero_gens, srels, gen_cols, pivots)
   # Compute map from elements of submodule into original module
   hmat = [nonzero_gens[gen_cols[i]].v[1, j] for i in 1:ngens(M) for j in 1:ngens(m)]
   f = ModuleHomomorphism(M, m, matrix(R, ngens(M), ngens(m), hmat))
   M.map = f
   return M, f
end

function Submodule(m::AbstractAlgebra.FPModule{T}, gens::Vector{Any}) where T <: RingElement
   length(gens) != 0 && error("Incompatible module elements")
   return Submodule(m, elem_type(m)[])
end
