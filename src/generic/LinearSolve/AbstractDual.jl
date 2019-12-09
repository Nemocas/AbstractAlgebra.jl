
# Unfortunately, the matrix interface has chained itself to commutative rings.
# It is going to be a royal pain to generalize the functionality to non-commutative
# rings. In the mean time, we will still support duality for the day that occurs.

abstract type OppositeRingFlag <: NCRing end
abstract type OppositeEltFlag <: NCRingElem end

# NOTE: Mutability is important for the constructors below.
mutable struct BasicOppositeEltFlag{T} <: NCRingElem
    elt::T
end

#@inline function BasicOppositeEltFlag(a::T) where T 
#    return BasicOppositeEltFlag{T}(a)
#end


@inline function element(a::BasicOppositeEltFlag)
    return a.elt
end

function set_element!(a::BasicOppositeEltFlag, b)
    a.elt = b
    return a
end

# Ensure `Dual(Dual(x)) === x`. This is why we need mutable struct.
@inline function (::Type{T})(Oelt::T) where T <: OppositeEltFlag
    return element(Oelt)
end

# Notice that the @inlines forces the compiler to compile Dual(Dual(a)) to just be a.

###############################################################################
#
#   Opposite Ring Elt interface
#
###############################################################################

function add!(c::OppositeEltFlag, a::OppositeEltFlag, b::OppositeEltFlag)
    celt = element(c)
    celt = add!(celt, element(a), element(b))
    c = set_element!(c,celt)
    return c
end

function mul!(c::OppositeEltFlag, a::OppositeEltFlag, b::OppositeEltFlag)
    celt = element(c)
    celt = add!(celt, element(b), element(a))
    c = set_element!(c,celt)
    return c
end

function +(a::OppositeEltFlag, b::OppositeEltFlag)
    c = deepcopy(a)
    return add!(c,a,b)
end

function *(a::OppositeEltFlag, b::OppositeEltFlag)
    c = deepcopy(a)
    return mul!(c,a,b)
end

function reduce!(a::OppositeEltFlag)
    elt = element(a)
    elt = reduce!(elt)
    return set_element!(a,elt)
end

function zero!(a::OppositeEltFlag)
    elt = element(a)
    elt = zero!(a)
    return set_element!(a,elt)
end

function divexact_right(a::OppositeEltFlag, b::OppositeEltFlag)
    elt_a = element(a)
    elt_b = element(b)
    celt = divexact_left(elt_b, elt_a)
    return typeof(a)(celt)
end

function divexact_left(a::OppositeEltFlag, b::OppositeEltFlag)
    elt_a = element(a)
    elt_b = element(b)
    celt = divexact_right(elt_b, elt_a)
    return typeof(a)(celt)
end

# Generally you get the point. For now, it is very easy to implement this interface to
# the extent we need it.

###############################################################################
#
#   Opposite Elt Flag (commutative rings)
#
###############################################################################

# Indeed, all our hard work above was for nothing.
function BasicOppositeEltFlag(a::RingElement)
    return a
end

###############################################################################
#
#   TransposeIndexDual
#
###############################################################################

mutable struct TransposeIndexDual{T} <: MatElem{T}
    elt::MatElem{T}
end

mutable struct MatrixOfOpposites{T} <: MatElem{T}
    elt::MatElem{T}
end


###############################################################################
#
#   Fancy silliness on top of the silliness.
#
###############################################################################

#@inline function TransposeIndexDual(elt::MatElem{T}) where T
#    return TransposeIndexDual{T}(elt)
#end

@inline function element(D::TransposeIndexDual)
    return D.elt
end

@inline function TransposeIndexDual(D::TransposeIndexDual)
    return element(D)
end

function set_element!(D, elt)
    D.elt = elt
    return D
end

####
# Implement a light version of the matrix interface.

INDEX_TYPES = Union{Colon, Int64, AbstractArray{Int64,1}}

function getindex(D::TransposeIndexDual, I::Int, J::Int)
    return BasicOppositeEltFlag(getindex(element(D), J, I))
end

function setindex!(D::TransposeIndexDual, a::OppositeEltFlag, I::INDEX_TYPES, J::INDEX_TYPES)
    elt = element(D)
    elt = setindex!(elt, J, I, typeof(a)(a))
    D = set_element!(D, elt)
    return D
end

# Catch the trivial case with commutative rings.
function setindex!(D::TransposeIndexDual, a::RingElement, I::INDEX_TYPES, J::INDEX_TYPES)
    elt = element(D)
    elt = setindex!(elt, a, J, I)
    return 
end

function ncols(D::TransposeIndexDual)
    return nrows(element(D))
end

function nrows(D::TransposeIndexDual)
    return ncols(element(D))
end

function size(D::TransposeIndexDual)
    a,b = size(element(D))
    return (b,a)
end

function base_ring(D::TransposeIndexDual)
    return base_ring(element(D))
end


function Base.view(D::TransposeIndexDual, I::Vararg{Any, N}) where N
    elt = element(D)
    elt_view = view(elt, I...)
    Dview = TransposeIndexDual(elt_view)
    return Dview
end

# Change the memory layout of the element of `D` so that the
# get/set index is accessed in column major order. As you can see,
# the correct solution is to take the dual of the dual of the dual.
function column_major_access_form(D::TransposeIndexDual)
    elt = element(D)
    new_elt = transpose(elt)
    return MatrixOfOpposites(new_elt)
end

function column_major_access_form(D::TransposeIndexDual{T}) where T <: RingElement
    elt = element(D)
    new_elt = transpose(elt)
    return new_elt
end

###############################################################################
#
#   Similar
#
###############################################################################


function _similar(x::TransposeIndexDual{T}, R::Ring, r::Int, c::Int) where T <: RingElement
    TT = elem_type(R)
    M = Matrix{TT}(undef, (c, r))
    z = x isa AbstractAlgebra.MatElem ? MatSpaceElem{TT}(M) : MatAlgElem{TT}(M)
    z.base_ring = R
   return TransposeIndexDual(z)
end


###############################################################################
#
#   More fancy silliness on top of the silliness.
#
###############################################################################

#@inline function MatrixOfOpposites(elt::MatElem{T}) where T
#    return MatrixOfOpposites{T}(elt)
#end

@inline function element(D::MatrixOfOpposites)
    return D.elt
end

@inline function MatrixOfOpposites(D::MatrixOfOpposites)
    return element(D)
end

function getindex(D::MatrixOfOpposites, I::INDEX_TYPES, J::INDEX_TYPES)
    return BasicOppositeEltFlag(getindex(element(D), J, I))
end

function setindex!(D::MatrixOfOpposites, I, J, a::OppositeEltFlag)
    elt = element(D)
    elt = setindex!(elt, I, J, element(a))
    D = set_element!(D, elt)
    return D
end

# Catch the trivial case with commutative rings.
function setindex!(D::MatrixOfOpposites, I::INDEX_TYPES, J::INDEX_TYPES, a::RingElement)
    elt = element(D)
    elt = setindex!(elt, I, J, a)
    return set_element!(D,elt)
end
