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

ngens(N::QuotientModule{T}) where T <: RingElement = length(N.gens)

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
   parent(v1) != parent(v2) && error("Incompatible module elements")
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, N::QuotientModule{T}) where T <: RingElement
   println(io, "Quotient module of:")
   show(io, N.m)
   println(io, "")
   println(io, " with relations:")
   show(io, N.rels)
end

function show(io::IO, N::QuotientModule{T}) where T <: FieldElement
   println(io, "Quotient space of:")
   show(io, N.m)
   println(io, "")
   println(io, " with relations:")
   show(io, N.rels)
end

function show(io::IO, v::quotient_module_elem)
   print(io, "(")
   len = ngens(parent(v))
   for i = 1:len - 1
      print(io, v.v[1, i])
      print(io, ", ")
   end
   if len > 0
      print(io, v.v[1, len])
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

function reduce_mod_rels(v::AbstractAlgebra.MatElem{T}, N::QuotientModule{T}) where T <: RingElement
   rels = N.rels # relation vectors (in terms of supermodule generators)
   gens = N.gens # columns of rels corresponding to gens of N
   culled = N.culled # list of indices of relations with non-unit pivot
   pivots = N.pivots # index of pivot in culled relations
   v = deepcopy(v) # don't destroy input
   geni = 1
   t = base_ring(N)()
   for i in 1:length(culled) # for each culled relation
      culli = culled[i] # index of culled relation in rels vector
      w = rels[culli].v # get next culled row as vector w
      col = pivots[i] # column of pivot in culled row
      # find index in gens of pivot column
      while gens[geni] < col
         geni += 1
      end
      d = w[1, col] # ring element in pivot column of culled row
      q, v[1, geni] = divrem(v[1, geni], d) # reduce entry in v by pivot
      # v = v - q*rels[culli], skipping zero cols
      q = -q
      for j = geni + 1:length(gens)
         t = mul!(t, q, w[1, gens[j]])
         v[1, j] = addeq!(v[1, j], t)
      end
   end
   return v 
end

function (N::QuotientModule{T})(v::Vector{T}) where T <: RingElement
   length(v) != ngens(N) && error("Length of vector does not match number of generators")
   mat = matrix(base_ring(N), 1, length(v), v)
   mat = reduce_mod_rels(mat, N)
   return quotient_module_elem{T}(N, mat)
end

function (N::QuotientModule{T})(v::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(v) != ngens(N) && error("Length of vector does not match number of generators")
   nrows(v) != 1 && ("Not a vector in quotient_module_elem constructor")
   v = reduce_mod_rels(v, N)
   return quotient_module_elem{T}(N, v)
end

###############################################################################
#
#   QuotientModule constructor
#
###############################################################################

@doc Markdown.doc"""
    QuotientModule(m::AbstractAlgebra.Module{T}, sub::Submodule{T}) where T <: RingElement
> Return the quotient of the module `m` by the module `sub`, which must have
> been constructed as a submodule of `m`.
"""
function QuotientModule(m::AbstractAlgebra.Module{T}, sub::Submodule{T}) where T <: RingElement
   supermodule(sub) !== m && error("Not a submodule in QuotientModule constructor") 
   M = QuotientModule{T}(m, sub.gens)
   G = gens(m)
   f = map_from_func(M, m, x -> sum(x.v[1, i]*G[M.gens[i]] for i in 1:ncols(x.v)))
   M.map = f
   return M, f
end

