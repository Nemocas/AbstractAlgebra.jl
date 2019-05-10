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

function check_parent(v1::quotient_module_elem{T}, v2::quotient_module_elem{T}) where T <: RingElement
   parent(v1) !== parent(v2) && error("Incompatible module elements")
end

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
   if length(relations(N)) == 0
      println("no relations")
   else
      println(io, "relations:")
      rels = [string(v) for v in relations(N)]
      print(IOContext(io, :compact => true), join(rels, ", "))
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
   N = parent(v)
   return N(v.v*c)
end

function *(v::quotient_module_elem{T}, c::U) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(v.v*c)
end

function *(c::T, v::quotient_module_elem{T}) where T <: RingElem
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

function ==(m::quotient_module_elem{T}, n::quotient_module_elem{T}) where T <: RingElement
   check_parent(m, n)
   return m.v == n.v
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function reduce_mod_rels(v::AbstractAlgebra.MatElem{T}, rels::Vector{<:AbstractAlgebra.MatElem{T}}) where T <: RingElement
   R = base_ring(v)
   v = deepcopy(v) # don't destroy input
   i = 1
   t1 = R()
   for rel in rels # for each relation
      while iszero(rel[1, i])
         i += 1
      end
      q, v[1, i] = AbstractAlgebra.divrem(v[1, i], rel[1, i])
      q = -q
      for j = i + 1:ncols(v)
         t1 = mul!(t1, q, rel[1, j])
         v[1, j] = addeq!(v[1, j], t1)
      end
      i += 1
   end
   return v 
end

function (N::QuotientModule{T})(v::Vector{T}) where T <: RingElement
   length(v) != ngens(N) && error("Length of vector does not match number of generators")
   mat = matrix(base_ring(N), 1, length(v), v)
   mat = reduce_mod_rels(mat, relations(N))
   return quotient_module_elem{T}(N, mat)
end

function (N::QuotientModule{T})(v::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(v) != ngens(N) && error("Length of vector does not match number of generators")
   nrows(v) != 1 && ("Not a vector in quotient_module_elem constructor")
   v = reduce_mod_rels(v, relations(N))
   return quotient_module_elem{T}(N, v)
end

###############################################################################
#
#   QuotientModule constructor
#
###############################################################################

function projection(v::AbstractAlgebra.MatElem{T}, rels::Vector{<:AbstractAlgebra.MatElem{T}}, N::QuotientModule{T}) where T <: RingElement
   R = base_ring(N)
   # reduce mod relations
   v = reduce_mod_rels(v, rels)
   # project down to quotient module
   r = zero_matrix(R, 1, ngens(N))
   for i = 1:ngens(N)
      r[1, i] = v[1, N.gen_cols[i]]
   end
   return quotient_module_elem{T}(N, r)
end

@doc Markdown.doc"""
    QuotientModule(m::AbstractAlgebra.FPModule{T}, sub::Submodule{T}) where T <: RingElement
> Return the quotient `M` of the module `m` by the module `sub` (which must
> have been constructed as a submodule of `m`) along with the canonical
> quotient map from `m` to `M`.
"""
function QuotientModule(m::AbstractAlgebra.FPModule{T}, sub::Submodule{T}) where T <: RingElement
   supermodule(sub) !== m && error("Not a submodule in QuotientModule constructor")
   nrels = ngens(sub)
   rels = Vector{dense_matrix_type(T)}(undef, nrels)
   gens = generators(sub)
   for i = 1:nrels
      rels[i] = gens[i].v
   end
   M = QuotientModule{T}(m, rels)
   f = map_from_func(m, M, x -> projection(x.v, rels, M))
   M.map = f
   return M, f
end

