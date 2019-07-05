###############################################################################
#
#   DirectSumModule.jl : Direct sums of modules
#
###############################################################################

export DirectSumModule, direct_sum_module_elem, summands

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{direct_sum_module_elem{T}}) where T <: RingElement = DirectSumModule{T}

elem_type(::Type{DirectSumModule{T}}) where T <: RingElement = direct_sum_module_elem{T}

parent(v::direct_sum_module_elem) = v.parent

base_ring(N::DirectSumModule{T}) where T <: RingElement = base_ring(N.m[1])

base_ring(v::direct_sum_module_elem{T}) where T <: RingElement = base_ring(v.parent)

ngens(N::DirectSumModule{T}) where T <: RingElement = sum(ngens(M) for M in N.m)

gens(N::DirectSumModule{T}) where T <: RingElement = [gen(N, i) for i = 1:ngens(N)]

function gen(N::DirectSumModule{T}, i::Int) where T <: RingElement
   R = base_ring(N)
   return N([(j == i ? one(R) : zero(R)) for j = 1:ngens(N)])
end

@doc Markdown.doc"""
    summands(M::DirectSumModule{T}) where T <: RingElement
> Return the modules that this module is a direct sum of.
"""
summands(M::DirectSumModule{T}) where T <: RingElement = M.m

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, N::DirectSumModule{T}) where T <: RingElement
   print(io, "DirectSumModule over ")
   print(IOContext(io, :compact => true), base_ring(N))
end

function show(io::IO, N::DirectSumModule{T}) where T <: FieldElement
   print(io, "Subspace over ")
   print(IOContext(io, :compact => true), base_ring(N))
end

function show(io::IO, v::direct_sum_module_elem)
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
#   Parent object call overload
#
###############################################################################

function (N::DirectSumModule{T})(v::Vector{T}) where T <: RingElement
   length(v) != ngens(N) && error("Length of vector does not match number of generators")
   mat = matrix(base_ring(N), 1, length(v), v)
   start = 1
   for i = 1:length(N.m)
      mat = reduce_mod_rels(mat, rels(N.m[i]), start)
      start += ngens(N.m[i])
   end
   return direct_sum_module_elem{T}(N, mat)
end

function (M::DirectSumModule{T})(a::Vector{Any}) where T <: RingElement
   length(a) != 0 && error("Incompatible element")
   return M(T[])
end

function (N::DirectSumModule{T})(v::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(v) != ngens(N) && error("Length of vector does not match number of generators")
   nrows(v) != 1 && ("Not a vector in direct_sum_module_elem constructor")
   start = 1
   for i = 1:length(N.m)
      v = reduce_mod_rels(v, rels(N.m[i]), start)
      start += ngens(N.m[i])
   end
   return direct_sum_module_elem{T}(N, v)
end

###############################################################################
#
#   DirectSum constructor
#
###############################################################################

function direct_sum_canonical_injection(m::AbstractAlgebra.FPModule{T}, D::DirectSumModule{T}, v::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
   R = base_ring(m)
   # Find starting point of the given module in the large vectors
   start = 0
   S = summands(D)
   i = 1
   while i < length(S) && S[i] !== m
         start += ngens(S[i])
         i += 1
   end
   # create embedded value
   newv = T[zero(R) for i in 1:ngens(D)]
   for i = 1:ngens(m)
      newv[i + start] = v[i]
   end
   matv = matrix(R, 1, length(newv), newv)
   return direct_sum_module_elem{T}(D, matv)
end

function direct_sum_canonical_projection(D::DirectSumModule{T}, m::U, v::AbstractAlgebra.FPModuleElem{T}) where {T <: RingElement, U <: AbstractAlgebra.FPModule{T}}
   # Find starting point of the given module in the large vectors
   R = base_ring(m)
   start = 0
   S = summands(D)
   i = 1
   while i < length(S) && S[i] !== m
         start += ngens(S[i])
         i += 1
   end
   # create projected value
   newv = T[v[i + start] for i in 1:ngens(m)]
   matv = matrix(R, 1, length(newv), newv)
   return elem_type(U)(m, matv)
end

@doc Markdown.doc"""
    DirectSum(m::Vector{AbstractAlgebra.FPModule{T}}) where T <: RingElement
> Return a tuple $M, f, g$ consisting of $M$ the direct sum of the modules `m`
> (supplied as a vector of modules), a vector $f$ of the canonical injections
> of the $m[i]$ into $M$ and a vector $g$ of the canonical projections from
> $M$ onto the $m[i]$.
"""
function DirectSum(m::Vector{<:AbstractAlgebra.FPModule{T}}) where T <: RingElement
   length(m) == 0 && error("Cannot take a direct sum of an empty vector of modules")
   # Check base rings are the same
   R = base_ring(m[1])
   for i = 2:length(m)
      base_ring(m[i]) != R && error("Incompatible modules")
   end
   # make vector of rels (only used by external interface, not internally)
   n = sum(ngens(m[i]) for i in 1:length(m))
   new_rels = Vector{dense_matrix_type(T)}(undef,
       sum(length(rels(m[i])) for i in 1:length(m)))
   rel = 1
   start = 0
   for i = 1:length(m)
      irels = rels(m[i])
      for j = 1:length(irels)
         new_rel = zero_matrix(R, 1, n)
         for k = 1:ncols(irels[j])
            new_rel[1, start + k] = irels[j][1, k]
         end
         new_rels[rel] = new_rel
         rel += 1
      end
   end
   # Construct DirectSumModule object
   M = DirectSumModule{T}(m, new_rels)
   # construct canonical injections and projections
   inj = Vector{ModuleHomomorphism{T}}(undef, length(m))
   pro = Vector{ModuleHomomorphism{T}}(undef, length(m))
   start = 0
   for i = 1:length(m)
      igens = ngens(m[i])
      mat1 = zero_matrix(R, igens, n)
      for j = 1:igens
         mat1[j, start + j] = one(R)
      end
      inj[i] = ModuleHomomorphism(m[i], M, mat1)
      mat2 = transpose(mat1)
      pro[i] = ModuleHomomorphism(M, m[i], mat2)
      # Override image_fns with fast versions that don't do matrix-vector mul
      inj[i].image_fn  = x -> direct_sum_canonical_injection(m[i], M, x)
      pro[i].image_fn = x -> direct_sum_canonical_projection(M, m[i], x)
      start += ngens(m[i])
   end
   M.inj = inj
   M.pro = pro
   return M, inj, pro
end

function DirectSum(vals::AbstractAlgebra.FPModule{T}...) where T <: RingElement
   return DirectSum([vals...])
end

