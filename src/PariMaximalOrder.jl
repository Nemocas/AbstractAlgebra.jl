###########################################################################################
#
#   PariMaximalOrders.jl : Pari maximal orders (using data from nfinit)
#
###########################################################################################

export PariMaximalOrder, approx, coprime_multiplier, intersect, bounded_ideals,
       numden, prime_decomposition, LLL_reduce, valuation, factor, factor_mul

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

type PariMaximalOrder{S, T} <: PariRing
   pari_nf::PariNumberField{S, T}
end

type PariIdeal{S, T} <: RingElem
   ideal::Ptr{Int}
   parent::PariMaximalOrder{S, T}

   function PariIdeal(a::Ptr{Int}, par::PariMaximalOrder)
      r = new(gclone(a), par)
      finalizer(r, _pari_ideal_clear_fn)
      return r
   end
end

_pari_ideal_clear_fn(a::PariIdeal) = gunclone(a.ideal)

###########################################################################################
#
#   Untyped low-level Pari functions
#
###########################################################################################

residue(data::Ptr{Int}) = pari_load(data, 3)

modulus(data::Ptr{Int}) = pari_load(data, 2)
   
function basistoalg(nf::Ptr{Int}, data::Ptr{Int})
   return ccall((:basistoalg, :libpari), Ptr{Int}, (Ptr{Int}, Ptr{Int}), nf, data)
end

function alg(nf::Ptr{Int}, data::Ptr{Int})
   s = basistoalg(nf, data)
   mods = residue(s)
   return mods
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function basis{S, T}(ord::PariMaximalOrder{S, T})
   data = pari_load(ord.pari_nf.data, 8)
   pol_type = PariPolyRing{PariRationalField, T}
   par = pol_type(PariQQ)
   return pari_vec{pol_type}(data, par)
end

###########################################################################################
#
#   Ideals of bounded norm
#
###########################################################################################

function bounded_ideals{S, T}(ord::PariMaximalOrder{S, T}, bound::Int)
   av = unsafe_load(avma, 1)
   vec = ccall((:ideallist0, :libpari), Ptr{Int},
               (Ptr{Int}, Int, Int), ord.pari_nf.data, bound, 4)
   vec_type = pari_vec{PariMaximalOrder{S, T}}
   A = Array(vec_type, bound)
   for i = 1:bound
      A[i] = vec_type(pari_load(vec, i + 1), ord)
   end
   unsafe_store!(avma, av, 1)
   return A
end

###########################################################################################
#
#   Prime decomposition
#
###########################################################################################

function prime_decomposition{S, T}(ord::PariMaximalOrder{S, T}, p::BigInt)
   av = unsafe_load(avma, 1)
   pr = pari(p)
   vec = ccall((:idealprimedec, :libpari), Ptr{Int},
               (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pr.d)
   return pari_vec{PariMaximalOrder{S, T}}(vec, ord)
end

function prime_decomposition{S, T}(ord::PariMaximalOrder{S, T}, p::Integer)
   return prime_decomposition(ord, BigInt(p))
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, ord::PariMaximalOrder)
   print(io, "Maximal order of Number Field over Rational Field with defining polynomial ")
   print(pol(ord.pari_nf), " with Z-basis ")
   pari_print(io, basis(ord).data)
end

function show(io::IO, id::PariIdeal)
   av = unsafe_load(avma, 1)
   id2 = ccall((:idealtwoelt, :libpari), Ptr{Int}, 
               (Ptr{Int}, Ptr{Int}), id.parent.pari_nf.data, id.ideal)
   pari_print(io, id2)
   unsafe_store!(avma, av, 1)
end

###########################################################################################
#
#   Ideal norm
#
###########################################################################################

function norm{S, T}(a::PariIdeal{S, T})
   av = unsafe_load(avma, 1)
   n = ccall((:idealnorm, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), a.parent.pari_nf.data, a.ideal)
   unsafe_store!(avma, av, 1)
   r = QQ()
   QQ!(r, n)
   return r
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

function *{S, T}(a::PariIdeal{S, T}, b::PariIdeal{S, T})
   av = unsafe_load(avma, 1)
   r = ccall((:idealmul, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), a.parent.pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(r, a.parent)
end

function +{S, T}(a::PariIdeal{S, T}, b::PariIdeal{S, T})
   av = unsafe_load(avma, 1)
   r = ccall((:idealadd, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), a.parent.pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(r, a.parent)
end

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{S, T}(a::PariIdeal{S, T}, n::Int)
   av = unsafe_load(avma, 1)
   r = ccall((:idealpows, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Int), a.parent.pari_nf.data, a.ideal, n)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(r, a.parent)
end

###########################################################################################
#
#   Ideal intersection
#
###########################################################################################

function intersect{S, T}(a::PariIdeal{S, T}, b::PariIdeal{S, T})
   av = unsafe_load(avma, 1)
   r = ccall((:idealintersect, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), a.parent.pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(r, a.parent)
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S, T}(a::PariIdeal{S, T}, b::PariIdeal{S, T})
   av = unsafe_load(avma, 1)
   r = ccall((:idealdiv, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), a.parent.pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(r, a.parent)
end

###########################################################################################
#
#   Inverse
#
###########################################################################################

function inv{S, T}(a::PariIdeal{S, T})
   av = unsafe_load(avma, 1)
   r = ccall((:idealinv, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), a.parent.pari_nf.data, a.ideal)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(r, a.parent)
end

###########################################################################################
#
#   LLL reduction
#
###########################################################################################

function LLL_reduce{S, T}(a::PariIdeal{S, T})
   av = unsafe_load(avma, 1)
   r = ccall((:idealred0, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Void}), a.parent.pari_nf.data, a.ideal, C_NULL)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(r, a.parent)
end

###########################################################################################
#
#   Bezout
#
###########################################################################################

function bezout{S, T}(a::PariIdeal{S, T}, b::PariIdeal{S, T})
   pari_nf = a.parent.pari_nf.data
   av = unsafe_load(avma, 1)
   st = ccall((:idealaddtoone, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), pari_nf, a.ideal, b.ideal)
   s = alg(pari_nf, pari_load(st, 2))
   t = alg(pari_nf, pari_load(st, 3))
   par = FmpqPolyRing{T}(QQ)
   pols = fmpq_poly!(par(), s)
   polt = fmpq_poly!(par(), t)
   unsafe_store!(avma, av, 1)
   return pols, polt
end

###########################################################################################
#
#   Numerator/denominator
#
###########################################################################################

function numden{S, T}(a::PariIdeal{S, T})
   pari_nf = a.parent.pari_nf
   av = unsafe_load(avma, 1)
   nd = ccall((:idealnumden, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), pari_nf.data, a.ideal)
   num = a.parent(pari_load(nd, 2))
   den = a.parent(pari_load(nd, 3))
   unsafe_store!(avma, av, 1)
   return num, den
end

###########################################################################################
#
#   Valuation
#
###########################################################################################

function valuation{S, T}(a::PariIdeal{S, T}, b::PariIdeal{S, T})
   pari_nf = a.parent.pari_nf
   av = unsafe_load(avma, 1)
   r = ccall((:idealval, :libpari), Int, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), pari_nf.data, a.ideal, b.ideal)
   unsafe_store!(avma, av, 1)
   return r
end

###########################################################################################
#
#   Ideal factorisation
#
###########################################################################################

function factor{S, T}(a::PariIdeal{S, T})
   av = unsafe_load(avma, 1)
   pari_fac = ccall((:idealfactor, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), a.parent.pari_nf.data, a.ideal)
   fac = PariFactor{PariMaximalOrder{S, T}}(pari_fac, a.parent)
   unsafe_store!(avma, av, 1)
   return fac
end

function factor_mul{S, T}(a::PariFactor{PariMaximalOrder{S, T}})
   av = unsafe_load(avma, 1)
   p = ccall((:idealfactorback, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Void}, Int), a.parent.pari_nf.data, a.data, C_NULL, 0)
   r = a.parent(p)
   unsafe_store!(avma, av, 1)
   return r
end

###########################################################################################
#
#   Approximation
#
###########################################################################################

function approx{S, T}(a::PariIdeal{S, T})
   pari_nf = a.parent.pari_nf.data
   par = a.parent
   av = unsafe_load(avma, 1)
   a = ccall((:idealappr, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}), pari_nf, a.ideal)
   r = alg(pari_nf, a)
   pol = FmpqPolyRing{T}(QQ)()
   fmpq_poly!(pol, r)
   unsafe_store!(avma, av, 1)
   return pol
end

###########################################################################################
#
#   Coprime ideal construction
#
###########################################################################################

function coprime_multiplier{S, T}(a::PariIdeal{S, T}, b::PariIdeal{S, T})
   pari_nf = a.parent.pari_nf.data
   par = a.parent
   av = unsafe_load(avma, 1)
   m = ccall((:idealcoprime, :libpari), Ptr{Int}, 
             (Ptr{Int}, Ptr{Int}, Ptr{Int}), pari_nf, a.ideal, b.ideal)
   r = alg(pari_nf, m)
   pol = FmpqPolyRing{T}(QQ)()
   fmpq_poly!(pol, r)
   unsafe_store!(avma, av, 1)
   return pol
end

###########################################################################################
#
#   Parent call object overloads
#
###########################################################################################

function Base.call{S, T}(ord::PariMaximalOrder{S, T}, ideal::Ptr{Int})
   return PariIdeal{S, T}(ideal, ord)
end

function Base.call{S, T}(ord::PariMaximalOrder{S, T}, b::fmpq_poly)
   av = unsafe_load(avma, 1)
   ideal = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(b).d)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(ideal, ord)
end

function Base.call{S, T}(ord::PariMaximalOrder{S, T}, b::nf_elem)
   av = unsafe_load(avma, 1)
   par = b.parent.pol.parent
   ideal = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(par(b)).d)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(ideal, ord)
end

function Base.call{S, T}(ord::PariMaximalOrder{S, T}, b::BigInt)
   av = unsafe_load(avma, 1)
   ideal = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(b).d)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(ideal, ord)
end

function Base.call{S, T}(ord::PariMaximalOrder{S, T}, b::Integer)
   av = unsafe_load(avma, 1)
   ideal = ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.pari_nf.data, pari(BigInt(b)).d)
   unsafe_store!(avma, av, 1)
   return PariIdeal{S, T}(ideal, ord)
end

###########################################################################################
#
#   MaximalOrder constructor
#
###########################################################################################

function MaximalOrder{S, T}(nf::NfNumberField{S, T})
   return PariMaximalOrder{S, T}(PariNumberField{S, T}(nf))
end
