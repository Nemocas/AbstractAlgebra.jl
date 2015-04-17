###############################################################################
#
#   PariMaximalOrders.jl : Pari maximal orders (using data from nfinit)
#
###############################################################################

export PariMaximalOrder

###############################################################################
#
#   Types and memory management
#
###############################################################################

type PariMaximalOrder <: PariRing
   pari_nf::PariNumberField
end

type PariMaximalOrderElem <: RingElem
   data::Ptr{Int}
   parent::PariMaximalOrder

   function PariMaximalOrderElem(a::Ptr{Int}, par::PariMaximalOrder)
      r = new(gclone(a), par)
      finalizer(r, _pari_maximal_order_elem_clear_fn)
      return r
   end
end

_pari_maximal_order_elem_clear_fn(a::PariMaximalOrderElem) = gunclone(a.data)

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
   pol_type = PariPolyRing{PariRationalField}
   par = pol_type(PariQQ, var(parent(ord.pari_nf.nf.pol)))
   return pari_vec{pol_type}(data, par)
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

function show(io::IO, el::PariMaximalOrderElem)
   av = unsafe_load(avma, 1)
   pari_print(io, el.data)
   unsafe_store!(avma, av, 1)
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function Base.call(ord::PariMaximalOrder, b::fmpq_poly)
   av = unsafe_load(avma, 1)
   data = ccall((:algtobasis, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(b).d)
   unsafe_store!(avma, av, 1)
   return PariMaximalOrderElem(data, ord)
end

function Base.call(ord::PariMaximalOrder, b::nf_elem)
   av = unsafe_load(avma, 1)
   par = b.parent.pol.parent
   data = ccall((:algtobasis, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(par(b)).d)
   unsafe_store!(avma, av, 1)
   return PariMaximalOrderElem(data, ord)
end

function Base.call(ord::PariMaximalOrder, b::fmpz)
   av = unsafe_load(avma, 1)
   data = ccall((:algtobasis, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(b).d)
   unsafe_store!(avma, av, 1)
   return PariMaximalOrderElem(data, ord)
end

function Base.call(ord::PariMaximalOrder, b::Integer)
   av = unsafe_load(avma, 1)
   data = ccall((:algtobasis, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(ZZ(b)).d)
   unsafe_store!(avma, av, 1)
   return PariMaximalOrderElem(data, ord)
end

###############################################################################
#
#   MaximalOrder constructor
#
###############################################################################

function MaximalOrder(nf::NfNumberField)
   return PariMaximalOrder(PariNumberField(nf))
end
