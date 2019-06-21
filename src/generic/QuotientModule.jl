###############################################################################
#
#   QuotientModule.jl : Quotients of modules by submodules
#
###############################################################################

export QuotientModule, quotient_module_elem

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{quotient_module_elem{T}}) where T <: RingElement = QuotientModule{T}

elem_type(::Type{QuotientModule{T}}) where T <: RingElement = quotient_module_elem{T}

parent(v::quotient_module_elem) = v.parent

base_ring(N::QuotientModule{T}) where T <: RingElement = N.base_ring

base_ring(v::quotient_module_elem{T}) where T <: RingElement = base_ring(v.parent)

ngens(N::QuotientModule{T}) where T <: RingElement = length(N.gen_cols)

gens(N::QuotientModule{T}) where T <: RingElement = [gen(N, i) for i = 1:ngens(N)]

function gen(N::QuotientModule{T}, i::Int) where T <: RingElement
   R = base_ring(N)
   mat = matrix(R, 1, ngens(N),
                [(j == i ? one(R) : zero(R)) for j = 1:ngens(N)])
   return quotient_module_elem{T}(N, mat)
end

@doc Markdown.doc"""
    supermodule(M::QuotientModule{T}) where T <: RingElement
> Return the module that this module is a quotient of.
"""
supermodule(M::QuotientModule{T}) where T <: RingElement = M.m

###############################################################################
#
#   String I/O
#
###############################################################################

function show_gens_rels(io::IO, N::AbstractAlgebra.FPModule{T}) where T <: RingElement
   print(io, " with ", ngens(N), " generator")
   if ngens(N) == 1
      print(io, " and ")
   else
      print(io, "s and ")
   end
   if length(rels(N)) == 0
      println(io, "no relations")
   else
      println(io, "relations:")
      Nrels = [string(v) for v in rels(N)]
      print(IOContext(io, :compact => true), join(Nrels, ", "))
   end
end

function show(io::IO, N::QuotientModule{T}) where T <: RingElement
   print(io, "Quotient module over ")
   print(IOContext(io, :compact => true), base_ring(N))
   show_gens_rels(io, N)
end

function show(io::IO, N::QuotientModule{T}) where T <: FieldElement
   println(io, "Quotient space over:")
   print(IOContext(io, :compact => true), base_ring(N))
   show_gens_rels(io, N)
end

function show(io::IO, v::quotient_module_elem)
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

function -(v::quotient_module_elem{T}) where T <: RingElement
   N = parent(v)
   return N(-v.v)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(v1::quotient_module_elem{T}, v2::quotient_module_elem{T}) where T <: RingElement
   check_parent(v1, v2)
   N = parent(v1)
   return N(v1.v + v2.v)
end

function -(v1::quotient_module_elem{T}, v2::quotient_module_elem{T}) where T <: RingElement
   check_parent(v1, v2)
   N = parent(v1)
   return N(v1.v - v2.v)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(v::quotient_module_elem{T}, c::T) where T <: RingElem
   base_ring(v) != parent(c) && error("Incompatible rings")
   N = parent(v)
   return N(v.v*c)
end

function *(v::quotient_module_elem{T}, c::U) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(v.v*c)
end

function *(c::T, v::quotient_module_elem{T}) where T <: RingElem
   base_ring(v) != parent(c) && error("Incompatible rings")
   N = parent(v)
   return N(c*v.v)
end

function *(c::U, v::quotient_module_elem{T}) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(c*v.v)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(M::QuotientModule{T}, N::QuotientModule{T}) where T <: RingElement
   check_parent(M, N)
   return M.m == N.m && M.rels == N.rels
end

function ==(m::quotient_module_elem{T}, n::quotient_module_elem{T}) where T <: RingElement
   check_parent(m, n)
   return m.v == n.v
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

# Reduce the vector v by the relations vrels starting at column start of v
function reduce_mod_rels(v::AbstractAlgebra.MatElem{T}, vrels::Vector{<:AbstractAlgebra.MatElem{T}}, start::Int) where T <: RingElement
   R = base_ring(v)
   v = deepcopy(v) # don't destroy input
   i = 1
   t1 = R()
   for k = 1:length(vrels) # for each relation
      rel = vrels[k]
      while iszero(rel[1, i])
         i += 1
      end
      q, v[1, start + i - 1] = AbstractAlgebra.divrem(v[1, start + i - 1], rel[1, i])
      q = -q
      for j = i + 1:ncols(rel)
         t1 = mul!(t1, q, rel[1, j])
         v[1, start + j - 1] = addeq!(v[1, start + j - 1], t1)
      end
      i += 1
   end
   return v 
end

function (N::QuotientModule{T})(v::Vector{T}) where T <: RingElement
   length(v) != ngens(N) && error("Length of vector does not match number of generators")
   mat = matrix(base_ring(N), 1, length(v), v)
   mat = reduce_mod_rels(mat, rels(N), 1)
   return quotient_module_elem{T}(N, mat)
end

function (M::QuotientModule{T})(a::Vector{Any}) where T <: Union{RingElement, NCRingElem}
   length(a) != 0 && error("Incompatible element")
   return M(T[])
end

function (N::QuotientModule{T})(v::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(v) != ngens(N) && error("Length of vector does not match number of generators")
   nrows(v) != 1 && ("Not a vector in quotient_module_elem constructor")
   v = reduce_mod_rels(v, rels(N), 1)
   return quotient_module_elem{T}(N, v)
end

###############################################################################
#
#   QuotientModule constructor
#
###############################################################################

function projection(v::AbstractAlgebra.MatElem{T}, crels::AbstractAlgebra.MatElem{T}, N::QuotientModule{T}) where {T <: RingElement}
   R = base_ring(N)
   # remove zero rows
   nr = nrows(crels)
   while nr > 0 && iszero_row(crels, nr)
      nr -= 1
   end
   # put into row vectors
   vrels = Vector{dense_matrix_type(T)}(undef, nr)
   for i = 1:nr
      vrels[i] = matrix(R, 1, ncols(crels), [crels[i, j] for j in 1:ncols(crels)])
   end
   # reduce mod relations
   v = reduce_mod_rels(v, vrels, 1)
   # project down to quotient module
   r = zero_matrix(R, 1, ngens(N))
   for i = 1:ngens(N)
      r[1, i] = v[1, N.gen_cols[i]]
   end
   return r
end

function compute_combined_rels(m::AbstractAlgebra.FPModule{T}, srels::Vector{S}) where {T <: RingElement, S <: AbstractAlgebra.MatElem{T}}
   # concatenate relations in m and new rels
   R = base_ring(m)
   old_rels = rels(m)
   combined_rels = zero_matrix(R, length(old_rels) + length(srels), ngens(m))
   for i = 1:length(old_rels)
      for j = 1:ngens(m)
         combined_rels[i, j] = old_rels[i][1, j]
      end
   end
   for i = 1:length(srels)
      for j = 1:ngens(m)
         combined_rels[i + length(old_rels), j] = srels[i][1, j]
      end
   end
   # compute the hnf/rref of the combined relations
   combined_rels = reduced_form(combined_rels)
   return combined_rels
end

function QuotientModule(m::AbstractAlgebra.FPModule{T}, sub::Submodule{T}) where T <: RingElement
   !issubmodule(m, sub) && error("Not a submodule in QuotientModule constructor")
   R = base_ring(m)
   if sub === m # quotient of submodule by itself
      srels = [v.v for v in gens(sub)]
      combined_rels = compute_combined_rels(m, srels)
      M = QuotientModule{T}(m, combined_rels)
      f = ModuleHomomorphism(m, M,
          matrix(R, ngens(m), 0, []))
   else
      G = generators(sub)
      S = sub
      if supermodule(S) !== m
         while supermodule(S) !== m
            G = elem_type(typeof(supermodule(S.m)))[S.m.map(v) for v in G]
            S = supermodule(S)
         end
         sub, v = Submodule(m, G)
         G = generators(sub)
      end
      nrels = ngens(sub)
      srels = Vector{dense_matrix_type(T)}(undef, nrels)
      for i = 1:nrels
         srels[i] = G[i].v
      end
      combined_rels = compute_combined_rels(m, srels)
      M = QuotientModule{T}(m, combined_rels)
      hvecs = [projection(x.v, combined_rels, M) for x in gens(m)]
      hmat = [hvecs[i][1, j] for i in 1:ngens(m) for j in 1:ngens(M)]
      f = ModuleHomomorphism(m, M, matrix(R, ngens(m), ngens(M), hmat))
   end
   M.map = f
   return M, f
end

@doc Markdown.doc"""
    QuotientModule(m::AbstractAlgebra.FPModule{T}, sub::AbstractAlgebra.FPModule{T}) where T <: RingElement
> Return the quotient `M` of the module `m` by the module `sub` (which must
> have been (transitively) constructed as a submodule of `m` or be `m` itself)
> along with the canonical quotient map from `m` to `M`.
"""
function QuotientModule(m::AbstractAlgebra.FPModule{T}, sub::AbstractAlgebra.FPModule{T}) where T <: RingElement
   # The only case we need to deal with here is where `m == sub`. In all other
   # cases, sub will be of type Submodule.
   m !== sub && error("Not a submodule in QuotientModule constructor")
   srels = [v.v for v in gens(sub)]
   combined_rels = compute_combined_rels(m, srels)
   M = QuotientModule{T}(m, combined_rels)
   f = ModuleHomomorphism(m, M, matrix(R, ngens(m), 0, []))
   M.map = f
   return M, f   
end
