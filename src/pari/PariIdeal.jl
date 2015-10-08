###############################################################################
#
#   PariIdeal.jl : Pari ideals in the maximal order (for now)
#
###############################################################################

export approx, coprime_multiplier, intersect, bounded_ideals, numden, 
       prime_decomposition, valuation, factor, factor_mul,
       PariIdealSet, PariIdeal, PariFactor

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent(a::PariIdeal) = a.parent

function check_parent(a::PariIdeal, b::PariIdeal) 
   a.parent != b.parent && error("Incompatible ideals in ideal operation")
end

###############################################################################
#
#   Ideals of bounded norm
#
###############################################################################

function bounded_ideals(ord::PariMaximalOrder, bound::Int)
   av = unsafe_load(avma, 1)
   vec = ccall((:ideallist0, :libpari), Ptr{Int},
               (Ptr{Int}, Int, Int), ord.pari_nf.data, bound, 4)
   vec_type = pari_vec{PariIdeal}
   A = Array(vec_type, bound)
   for i = 1:bound
      A[i] = vec_type(pari_load(vec, i + 1), PariIdealSet(ord))
   end
   unsafe_store!(avma, av, 1)
   return A
end

###############################################################################
#
#   Prime decomposition
#
###############################################################################

function prime_decomposition(ord::PariMaximalOrder, p::fmpz)
   av = unsafe_load(avma, 1)
   pr = pari(p)
   vec = ccall((:idealprimedec, :libpari), Ptr{Int},
               (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pr.d)
   return pari_vec{PariIdeal}(vec, PariIdealSet(ord))
end

function prime_decomposition(ord::PariMaximalOrder, p::Integer)
   return prime_decomposition(ord, fmpz(p))
end

###############################################################################
#
#   AbstractString{} I/O
#
###############################################################################

function show(io::IO, s::PariIdealSet)
   print(io, "Set of ideals of ")
   print(io, s.order)
end

function show(io::IO, id::PariIdeal)
   av = unsafe_load(avma, 1)
   id2 = ccall((:idealtwoelt, :libpari), Ptr{Int}, 
               (Ptr{Int}, Ptr{Int}), id.parent.order.pari_nf.data, id.ideal)
   pari_print(io, id2)
   unsafe_store!(avma, av, 1)
end

###############################################################################
#
#   Ideal norm
#
###############################################################################

function norm(a::PariIdeal)
   av = unsafe_load(avma, 1)
   n = ccall((:idealnorm, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), a.parent.order.pari_nf.data, a.ideal)
   unsafe_store!(avma, av, 1)
   r = fmpq()
   fmpq!(r, n)
   return r
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function *(a::PariIdeal, b::PariIdeal)
   check_parent(a, b)
   av = unsafe_load(avma, 1)
   r = ccall((:idealmul, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), 
                   a.parent.order.pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal(r, a.parent)
end

function +(a::PariIdeal, b::PariIdeal)
   check_parent(a, b)
   av = unsafe_load(avma, 1)
   r = ccall((:idealadd, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), 
                   a.parent.order.pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal(r, a.parent)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::PariIdeal, b::PariIdeal)
   check_parent(a, b)
   av = unsafe_load(avma, 1)
   r1 = ccall((:idealhnf, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), 
                   a.parent.order.pari_nf.data, a.ideal)
   r2 = ccall((:idealhnf, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), 
                   a.parent.order.pari_nf.data, b.ideal)
   r = Bool(ccall((:gequal, :libpari), Int, 
                  (Ptr{Int}, Ptr{Int}), r1, r2))
   unsafe_store!(avma, av, 1)
   return r
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::PariIdeal, n::Int)
   av = unsafe_load(avma, 1)
   r = ccall((:idealpows, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Int), 
                   a.parent.order.pari_nf.data, a.ideal, n)
   unsafe_store!(avma, av, 1)
   return PariIdeal(r, a.parent)
end

###############################################################################
#
#   Ideal intersection
#
###############################################################################

function intersect(a::PariIdeal, b::PariIdeal)
   check_parent(a, b)
   av = unsafe_load(avma, 1)
   r = ccall((:idealintersect, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), 
                   a.parent.order.pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal(r, a.parent)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::PariIdeal, b::PariIdeal)
   check_parent(a, b)
   av = unsafe_load(avma, 1)
   r = ccall((:idealdiv, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), 
                a.parent.order.pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal(r, a.parent)
end

###############################################################################
#
#   Inverse
#
###############################################################################

function inv(a::PariIdeal)
   av = unsafe_load(avma, 1)
   r = ccall((:idealinv, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), a.parent.order.pari_nf.data, a.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal(r, a.parent)
end

###############################################################################
#
#   GCDX
#
###############################################################################

function gcdx(a::PariIdeal, b::PariIdeal)
   check_parent(a, b)
   pari_nf = a.parent.order.pari_nf.data
   av = unsafe_load(avma, 1)
   st = ccall((:idealaddtoone, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), pari_nf, a.ideal, b.ideal)
   s = alg(pari_nf, pari_load(st, 2))
   t = alg(pari_nf, pari_load(st, 3))
   par = FmpqPolyRing(var(a.parent.order.pari_nf.nf.pol.parent))
   pols = fmpq_poly!(par(), s)
   polt = fmpq_poly!(par(), t)
   unsafe_store!(avma, av, 1)
   return pols, polt
end

###############################################################################
#
#   Numerator/denominator
#
###############################################################################

function numden(a::PariIdeal)
   pari_nf = a.parent.order.pari_nf
   av = unsafe_load(avma, 1)
   nd = ccall((:idealnumden, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), pari_nf.data, a.ideal)
   num = a.parent(pari_load(nd, 2))
   den = a.parent(pari_load(nd, 3))
   unsafe_store!(avma, av, 1)
   return num, den
end

###############################################################################
#
#   Valuation
#
###############################################################################

function valuation(a::PariIdeal, b::PariIdeal)
   check_parent(a, b)
   pari_nf = a.parent.order.pari_nf
   av = unsafe_load(avma, 1)
   r = ccall((:idealval, :libpari), Int, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return r
end

###############################################################################
#
#   Ideal factorisation
#
###############################################################################

function factor(a::PariIdeal)
   av = unsafe_load(avma, 1)
   pari_fac = ccall((:idealfactor, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), a.parent.order.pari_nf.data, a.ideal)
   fac = PariFactor{pari_maximal_order_elem}(pari_fac, a.parent)
   unsafe_store!(avma, av, 1)
   return fac
end

function factor_mul(a::PariFactor{pari_maximal_order_elem})
   av = unsafe_load(avma, 1)
   p = ccall((:idealfactorback, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Void}, Int), 
                    a.parent.order.pari_nf.data, a.data, C_NULL, 0)
   r = a.parent(p)
   unsafe_store!(avma, av, 1)
   return r
end

###############################################################################
#
#   Approximation
#
###############################################################################

function approx(a::PariIdeal)
   pari_nf = a.parent.order.pari_nf.data
   par = a.parent
   av = unsafe_load(avma, 1)
   a1 = ccall((:idealappr, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), pari_nf, a.ideal)
   r = alg(pari_nf, a1)
   pol = FmpqPolyRing(FlintQQ, var(a.parent.order.pari_nf.nf.pol.parent))()
   fmpq_poly!(pol, r)
   unsafe_store!(avma, av, 1)
   return pol
end

###############################################################################
#
#   Coprime ideal construction
#
###############################################################################

function coprime_multiplier(a::PariIdeal, b::PariIdeal)
   check_parent(a, b)
   pari_nf = a.parent.order.pari_nf.data
   par = a.parent
   av = unsafe_load(avma, 1)
   m = ccall((:idealcoprime, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), pari_nf, a.ideal, b.ideal)
   r = alg(pari_nf, m)
   pol = FmpqPolyRing(FlintQQ, var(a.parent.order.pari_nf.nf.pol.parent))()
   fmpq_poly!(pol, r)
   unsafe_store!(avma, av, 1)
   return pol
end

###############################################################################
#
#   Ideal creation functions
#
###############################################################################

function PariIdeal(ord::PariMaximalOrder, b::fmpq_poly, args::fmpq_poly...)
   ord.pari_nf.nf.pol.parent != parent(b) && error("Incompatible maximal order and polynomial")
   av = unsafe_load(avma, 1)
   id1 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(b).d)
   for pol in args
      id2 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(pol).d)
      id1 = ccall((:idealadd, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), ord.pari_nf.data, id1, id2)
   end
   unsafe_store!(avma, av, 1)
   return PariIdeal(id1, PariIdealSet(ord))
end

function PariIdeal(ord::PariMaximalOrder, b::pari_maximal_order_elem, args::pari_maximal_order_elem...)
   parent(b) != ord && error("Unable to coerce maximal order element")
   av = unsafe_load(avma, 1)
   id1 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, b.data)
   for el in args
      parent(el) != ord && error("Unable to coerce maximal order element")
      id2 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, el.data)
      id1 = ccall((:idealadd, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), ord.pari_nf.data, id1, id2)
   end
   unsafe_store!(avma, av, 1)
   return PariIdeal(id1, PariIdealSet(ord))
end

function PariIdeal(ord::PariMaximalOrder, b::nf_elem, args::nf_elem...)
   ord.pari_nf.nf != parent(b) && error("Incompatible maximal order and number field element")
   av = unsafe_load(avma, 1)
   par = b.parent.pol.parent
   id1 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(par(b)).d)
   for el in args
      id2 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(par(el)).d)
      id1 = ccall((:idealadd, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), ord.pari_nf.data, id1, id2)
   end
   unsafe_store!(avma, av, 1)
   return PariIdeal(id1, PariIdealSet(ord))
end

function PariIdeal(ord::PariMaximalOrder, b::fmpz, args::fmpz...)
   av = unsafe_load(avma, 1)
   id1 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(b).d)
   for c in args
      id2 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(c).d)
      id1 = ccall((:idealadd, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), ord.pari_nf.data, id1, id2)
   end
   unsafe_store!(avma, av, 1)
   return PariIdeal(id1, PariIdealSet(ord))
end

function PariIdeal(ord::PariMaximalOrder, b::Integer, args::Integer...)
   av = unsafe_load(avma, 1)
   id1 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(fmpz(b)).d)
   for c in args
      id2 = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(fmpz(c)).d)
      id1 = ccall((:idealadd, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), ord.pari_nf.data, id1, id2)
   end
   unsafe_store!(avma, av, 1)
   return PariIdeal(id1, PariIdealSet(ord))
end

###############################################################################
#
#   Parent object overloads
#
###############################################################################

function Base.call(ord::PariIdealSet, id::Ptr{Int})
   return PariIdeal(id, PariIdealSet(ord))
end

