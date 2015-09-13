###############################################################################
#
#   PariMaximalOrders.jl : Pari maximal orders (using data from nfinit)
#
###############################################################################

export PariMaximalOrder, pari_maximal_order_elem, basistoalg, alg, basis

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent(a::pari_maximal_order_elem) = a.parent

###############################################################################
#
#   Untyped low-level Pari functions
#
###############################################################################

residue(data::Ptr{Int}) = pari_load(data, 3)

modulus(data::Ptr{Int}) = pari_load(data, 2)
   
function basistoalg(nf::Ptr{Int}, data::Ptr{Int})
   return ccall((:basistoalg, :libpari), Ptr{Int}, 
                (Ptr{Int}, Ptr{Int}), nf, data)
end

function alg(nf::Ptr{Int}, data::Ptr{Int})
   s = basistoalg(nf, data)
   mods = residue(s)
   return mods
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function basis(ord::PariMaximalOrder)
   data = pari_load(ord.pari_nf.data, 8)
   poly_type = pari_poly{pari_rat}
   polyring = PariPolyRing{pari_rat}
   par = polyring(PariQQ, var(parent(ord.pari_nf.nf.pol)))
   return pari_vec{poly_type}(data, par)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, ord::PariMaximalOrder)
   print(io, "Maximal order of Number Field over Rational Field")
   print(io, " with defining polynomial ")
   print(pol(ord.pari_nf), " with Z-basis ")
   pari_print(io, basis(ord).data)
end

function show(io::IO, el::pari_maximal_order_elem)
   av = unsafe_load(avma, 1)
   pari_print(io, el.data)
   unsafe_store!(avma, av, 1)
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function Base.call(ord::PariMaximalOrder)
   av = unsafe_load(avma, 1)
   data = ccall((:algtobasis, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(fmpz()).d)
   unsafe_store!(avma, av, 1)
   return pari_maximal_order_elem(data, ord)
end

function Base.call(ord::PariMaximalOrder, b::fmpq_poly)
   av = unsafe_load(avma, 1)
   data = ccall((:algtobasis, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(b).d)
   unsafe_store!(avma, av, 1)
   return pari_maximal_order_elem(data, ord)
end

function Base.call(ord::PariMaximalOrder, b::nf_elem)
   av = unsafe_load(avma, 1)
   par = b.parent.pol.parent
   data = ccall((:algtobasis, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(par(b)).d)
   unsafe_store!(avma, av, 1)
   return pari_maximal_order_elem(data, ord)
end

function Base.call(ord::PariMaximalOrder, b::fmpz)
   av = unsafe_load(avma, 1)
   data = ccall((:algtobasis, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(b).d)
   unsafe_store!(avma, av, 1)
   return pari_maximal_order_elem(data, ord)
end

function Base.call(ord::PariMaximalOrder, b::Integer)
   av = unsafe_load(avma, 1)
   data = ccall((:algtobasis, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(fmpz(b)).d)
   unsafe_store!(avma, av, 1)
   return pari_maximal_order_elem(data, ord)
end

function Base.call(ord::PariMaximalOrder, b::pari_maximal_order_elem)
   parent(b) != ord && error("Unable to coerce maximal order element")
   return b
end

###############################################################################
#
#   PariMaximalOrder constructor
#
###############################################################################

function PariMaximalOrder(nf::AnticNumberField)
   return PariMaximalOrder(PariNumberField(nf))
end
