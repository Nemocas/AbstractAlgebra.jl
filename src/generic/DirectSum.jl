###############################################################################
#
#   DirectSumModule.jl : Generic direct sums of modules
#
###############################################################################

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{DirectSumModuleElem{T}}) where T <: RingElement = DirectSumModule{T}

elem_type(::Type{DirectSumModule{T}}) where T <: RingElement = DirectSumModuleElem{T}

parent(v::DirectSumModuleElem) = v.parent

base_ring_type(::Type{DirectSumModule{T}}) where T <: RingElement = parent_type(T)

base_ring(N::DirectSumModule{T}) where T <: RingElement = base_ring(N.m[1])::base_ring_type(N)

dim(M::DirectSumModule{<:FieldElem}) = sum(dim(x) for x = M.m)

@attr Int number_of_generators(N::DirectSumModule{T}) where T <: RingElement = sum(ngens(M) for M in N.m)

gens(N::DirectSumModule{T}) where T <: RingElement = [gen(N, i) for i = 1:ngens(N)]

rank(M::DirectSumModule{T}) where T = sum(rank, summands(M); init=0)
vector_space_dim(M::DirectSumModule{<:FieldElem}) = rank(M)

function gen(N::DirectSumModule{T}, i::Int) where T <: RingElement
   @boundscheck 1 <= i <= ngens(N) || throw(ArgumentError("generator index is out of range"))
   R = base_ring(N)
   return N([(j == i ? one(R) : zero(R)) for j = 1:ngens(N)])
end

@doc raw"""
    summands(M::DirectSumModule{T}) where T <: RingElement

Return the modules that this module is a direct sum of.
"""
summands(M::DirectSumModule{T}) where T <: RingElement = M.m

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, N::DirectSumModule{T}) where T <: RingElement
   @show_name(io, N)
   @show_special(io, N)
   if is_terse(io)
     io = pretty(io)
     print(io, LowercaseOff(), "DirectSumModule")
   else
     io = pretty(io)
     print(io, LowercaseOff(), "DirectSumModule over ")
     print(terse(io), Lowercase(), base_ring(N))
   end
end

function show(io::IO, v::DirectSumModuleElem)
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
   return DirectSumModuleElem{T}(N, mat)
end

function (M::DirectSumModule{T})(a::Vector{Any}) where T <: RingElement
   length(a) != 0 && error("Incompatible element")
   return M(T[])
end

function (N::DirectSumModule{T})(v::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(v) != ngens(N) && error("Length of vector does not match number of generators")
   nrows(v) != 1 && ("Not a vector in DirectSumModuleElem constructor")
   start = 1
   for i = 1:length(N.m)
     if length(rels(N.m[i])) > 0
        v = reduce_mod_rels(v, rels(N.m[i]), start)
      end
      start += ngens(N.m[i])
   end
   return DirectSumModuleElem{T}(N, v)
end

function (M::DirectSumModule{T})(a::SubmoduleElem{T}) where T <: RingElement
   R = parent(a)
   base_ring(R) != base_ring(M) && error("Incompatible modules")
   return M(R.map(a))
end

function (M::DirectSumModule{T})(a::DirectSumModuleElem{T}) where T <: RingElement
   R = parent(a)
   R != M && error("Incompatible modules")
   return a
end

function (M::DirectSumModule{T})(a::Vector{U}) where U <: AbstractAlgebra.FPModuleElem{T} where T <: RingElement
   S = summands(M)
   if length(a) != length(S)
     error("need one entry per summand")
   end
   if any(x->parent(a[x]) != S[x], 1:length(a))
     error("Incompatible modules")
   end
   s = M.inj[1](a[1])
   for i = 2:length(a)
     s += M.inj[i](a[i])
   end
   return s
end

# Fallback for all other kinds of modules
function (M::DirectSumModule{T})(a::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
   error("Unable to coerce into given module")
end

###############################################################################
#
#   DirectSum constructor
#
###############################################################################

function direct_sum_injection(i::Int, D::DirectSumModule{T}, v::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
   S = summands(D)
   m = S[i]
   R = base_ring(m)
   # Find starting point of the given module in the large vectors
   start = sum(map(x->ngens(x)::Int, S[1:i-1]))
   # create embedded value
   X = zero(D)
   for j = 1:ngens(m)
      X.v[j + start] = v[j]
   end
   return X
end

function add_direct_sum_injection!(X::DirectSumModuleElem{T}, i::Int, v::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
   D = parent(X)
   S = summands(D)
   m = S[i]
   R = base_ring(m)
   # Find starting point of the given module in the large vectors
   start = sum(map(x->ngens(x)::Int, S[1:i-1]))
   # create embedded value
   for j = 1:ngens(m)
     X.v[j+start] += v[j]
   end
   return X
end


function AbstractAlgebra.canonical_injection(A::DirectSumModule, i::Int)
  B = summands(A)[i]
  return ModuleHomomorphism(B, A, inj_proj_mat(base_ring(A), ngens(B), ngens(A), sum(ngens(x) for x = summands(A)[1:i-1]; init = 0)+1))
end

AbstractAlgebra._number_of_direct_product_factors(A::DirectSumModule) = length(summands(A))

function direct_sum_projection(D::DirectSumModule{T}, i::Int, v::AbstractAlgebra.FPModuleElem{T}) where {T <: RingElement}
   # Find starting point of the given module in the large vectors
   S = summands(D)
   m = S[i]
   R = base_ring(m)
   start = 0
   degs = get_attribute(D, :degs)::Vector{Int}
   for j=1:i-1
     start += degs[j]
   end
#   start = sum(map(x->ngens(x)::Int, S[1:i-1]))
   # create projected value
   X = zero(m)
   for j=1:ngens(m)
     X.v[j] = v[j+start]
   end
   return X
   newv = T[v[j + start] for j in 1:ngens(m)]
   matv = matrix(R, 1, length(newv), newv)
   return elem_type(m)(m, matv)
end

function AbstractAlgebra.canonical_projection(A::DirectSumModule, i::Int)
  B = summands(A)[i]
  return ModuleHomomorphism(A, B, inj_proj_mat(base_ring(A), ngens(A), ngens(B), sum(ngens(x) for x = summands(A)[1:i-1]; init = 0)+1))
end

function direct_sum(m::Vector{<:AbstractAlgebra.FPModule{T}}) where T <: RingElement
   length(m) == 0 && error("Cannot take a direct sum of an empty vector of modules")
   # Check base rings are the same
   R = base_ring(m[1])
   for i = 2:length(m)
      base_ring(m[i]) != R && error("Incompatible modules")
   end
   # make vector of rels (only used by external interface, not internally)
   degs = map(ngens, m)
   n = sum(degs)
   new_rels = Vector{dense_matrix_type(T)}(undef,
       sum(length(rels(m[i])) for i in 1:length(m)))
   if length(new_rels) > 0
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
   end
   # Construct DirectSumModule object
   M = DirectSumModule{T}(m, new_rels)
   set_attribute!(M, :degs => degs)
   # construct injections and projections
   inj = Vector{ModuleHomomorphism{T}}(undef, length(m))
   pro = Vector{ModuleHomomorphism{T}}(undef, length(m))
   start = 0
   for i = 1:length(m)
      igens = degs[i]
      inj[i] = ModuleHomomorphism(m[i], M, inj_proj_mat(R, igens, n, start+1))
      pro[i] = ModuleHomomorphism(M, m[i], inj_proj_mat(R, n, igens, start+1))
      # Override image_fns with fast versions that don't do matrix-vector mul
      inj[i].image_fn  = x -> direct_sum_injection(i, M, x)
      pro[i].image_fn = x -> direct_sum_projection(M, i, x)
      start += igens
   end
   M.inj = inj
   M.pro = pro
   return M, inj, pro
end

function ModuleHomomorphism(D::DirectSumModule{T}, A::AbstractAlgebra.FPModule{T}, m::Vector{<:ModuleHomomorphism{T}}) where T <: RingElement
   S = summands(D)
   length(S) == length(m) || error("map array has wrong length")
   all(i->domain(m[i]) == S[i] && codomain(m[i]) == A, 1:length(S)) ||
                            error("modules and maps are not compatible")
   return ModuleHomomorphism(D, A, vcat([x.matrix for x = m]...))
end

function ModuleHomomorphism(D::DirectSumModule{T}, A::DirectSumModule{T}, m::Matrix{<:Any}) where T <: RingElement
   SD = summands(D)
   SA = summands(A)
   size(m) == (length(SD), length(SA)) || error("dimensions do not match")
   for i = 1:length(SD)
      for j = 1:length(SA)
         if m[i, j] == 0
            m[i, j] = ModuleHomomorphism(SD[i], SA[j], [zero(SA[j]) for x = 1:ngens(SD[i])])
         else
            isa(m[i, j], ModuleHomomorphism{T}) || error("matrix must contain only 0 and homs")
            domain(m[i, j]) === SD[i] && codomain(m[i, j]) === SA[j] ||
                                    error("modules and maps are not compatible")
         end
      end
   end

   return ModuleHomomorphism(D, A, transpose(hvcat(Tuple([length(SD) for i = 1:length(SA)]), map(x->transpose(x.matrix), m)...)))
end

function hom_direct_sum(M::DirectSumModule{T}, N::DirectSumModule{T}, mp::Vector{ModuleHomomorphism{T}}) where T
  @assert length(M.m) == length(mp) == length(N.m)
  return hom(M, N, cat(map(matrix, mp)..., dims = (1,2)))
end

function hom_direct_sum(A::DirectSumModule{T}, B::DirectSumModule{T}, M::Matrix{<:Map{<:AbstractAlgebra.FPModule{T}, <:AbstractAlgebra.FPModule{T}}}) where {T}
  pro = canonical_projections(A)
  im = canonical_injections(B)
  s = hom(A, B, [zero(B) for i = 1:dim(A)])
  for i=1:length(pro)
    for j=1:length(im)
      s += pro[i]*M[i,j]*im[j]
    end
  end
  return s
end

