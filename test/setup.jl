using AbstractAlgebra

using SparseArrays, LinearAlgebra
using AbstractAlgebra: mul! # disambiguate from LinearAlgebra.mul!

if VERSION < v"0.7.0-DEV.2004"
   using Base.Test
else
   using Test
end

using Random: Random, MersenneTwister, randsubseq, shuffle

const SETUP = true

const rng = MersenneTwister()

const RINGS = Dict(
   "exact ring"          => (ZZ,                 (-1000:1000,)),
   "exact field"         => (GF(7),              ()),
   "inexact ring"        => (RealField["t"][1],  (0:200, -1000:1000)),
   "inexact field"       => (RealField,          (-1000:1000,)),
   "non-integral domain" => (ResidueRing(ZZ, 6), (0:5,)),
   "fraction field"      => (QQ,                 (-1000:1000,)),
)

# Perm-test.jl
const IntTypes = [Int8, Int16, Int32, Int, UInt8, UInt16, UInt32, UInt, BigInt]

# MatrixAlgebra-test.jl, Matrix-test.jl
primes100 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
499, 503, 509, 521, 523, 541]

function randprime(n::Int)
   if n > 100 || n < 1
      throw(DomainError())
   end
   return primes100[rand(1:n)]
end

function istriu(A::Union{Generic.Mat, Generic.MatAlgElem})
   m = nrows(A)
   n = ncols(A)
   d = 0
   for c = 1:n
      for r = m:-1:1
         if !iszero(A[r,c])
            if r < d
               return false
            end
            d = r
            break
         end
      end
   end
   return true
end

function is_snf(A::Union{Generic.Mat, Generic.MatAlgElem})
   m = nrows(A)
   n = ncols(A)
   a = A[1,1]
   for i = 2:min(m,n)
      q, r = divrem(A[i,i], a)
      if !iszero(r)
         return false
      end
      a = A[i,i]
   end
   for i = 1:n
      for j = 1:m
         if i == j
            continue
         end
         if !iszero(A[j,i])
            return false
         end
      end
   end
   return true
end

function is_weak_popov(P::Union{Generic.Mat, Generic.MatAlgElem}, rank::Int)
   zero_rows = 0
   pivots = zeros(ncols(P))
   for r = 1:nrows(P)
      p = AbstractAlgebra.find_pivot_popov(P, r)
      if P[r,p] == 0
         zero_rows += 1
         continue
      end
      if pivots[p] != 0
         return false
      end
      pivots[p] = r
   end
   if zero_rows != nrows(P)-rank
      return false
   end
   return true
end

# Simulate user matrix type belonging to AbstractArray
# with getindex but no setindex!
struct MyTestMatrix{T} <: AbstractArray{T, 2}
   d::T
   dim::Int
end

Base.getindex(a::MyTestMatrix{T}, r::Int, c::Int) where T = a.d

Base.size(a::MyTestMatrix{T}) where T = a.dim, a.dim

# Simulate user Field, together with a specialized matrix type
# (like fmpz / fmpz_mat)
struct F2 <: AbstractAlgebra.Field end

Base.zero(::F2) = F2Elem(false)
Base.one(::F2) = F2Elem(true)
(::F2)() = F2Elem(false)

struct F2Elem <: AbstractAlgebra.FieldElem
   x::Bool
end

(::F2)(x::F2Elem) = x
Base.:-(x::F2Elem) = x
Base.:+(x::F2Elem, y::F2Elem) = F2Elem(x.x ⊻ y.x)
Base.inv(x::F2Elem) = x.x ? x : throw(DivideError())
Base.:*(x::F2Elem, y::F2Elem) = F2Elem(x.x * y.x)

Base.convert(::Type{F2Elem}, x::Integer) = F2Elem(x % Bool)
Base.:(==)(x::F2Elem, y::F2Elem) = x.x == y.x

AbstractAlgebra.parent_type(::Type{F2Elem}) = F2
AbstractAlgebra.elem_type(::Type{F2}) = F2Elem
AbstractAlgebra.parent(x::F2Elem) = F2()
AbstractAlgebra.mul!(x::F2Elem, y::F2Elem, z::F2Elem) = y * z
AbstractAlgebra.addeq!(x::F2Elem, y::F2Elem) = x + y
AbstractAlgebra.divexact(x::F2Elem, y::F2Elem) = y.x ? x : throw(DivideError())

struct F2Matrix <: AbstractAlgebra.MatElem{F2Elem}
   m::Generic.MatSpaceElem{F2Elem}
end

AbstractAlgebra.nrows(a::F2Matrix) = nrows(a.m)
AbstractAlgebra.ncols(a::F2Matrix) = ncols(a.m)
AbstractAlgebra.base_ring(::F2Matrix) = F2()

Base.getindex(a::F2Matrix, r::Int64, c::Int64) = a.m[r, c]
Base.setindex!(a::F2Matrix, x::F2Elem, r::Int64, c::Int64) = a.m[r, c] = x
Base.similar(x::F2Matrix, R::F2, r::Int, c::Int) = F2Matrix(similar(x.m, r, c))

function AbstractAlgebra.zero_matrix(R::F2, r::Int, c::Int)
   mat = Array{F2Elem}(undef, r, c)
   for i=1:r, j=1:c
      mat[i, j] = zero(R)
   end
   z = Generic.MatSpaceElem{F2Elem}(mat)
   z.base_ring = R
   return F2Matrix(z)
end

function AbstractAlgebra.matrix(R::F2, mat::AbstractMatrix{F2Elem})
   mat = convert(Matrix, mat)
   z = Generic.MatSpaceElem{F2Elem}(mat)
   z.base_ring = R
   return F2Matrix(z)
end

function AbstractAlgebra.matrix(R::F2, r::Int, c::Int, mat::AbstractMatrix{F2Elem})
   AbstractAlgebra._check_dim(r, c, mat)
   matrix(R, mat)
end

# Modules-test.jl

function rand_module(R::AbstractAlgebra.Ring, vals...)
   rk = rand(0:5)
   M = FreeModule(R, rk)
   levels = rand(0:3)
   for i = 1:levels
      if ngens(M) == 0
         break
      end
      G = [rand(M, vals...) for i in 1:rand(1:ngens(M))]
      S, f = sub(M, G)
      if rand(1:2) == 1
         M, f = quo(M, S)
      else
         M = S
      end
   end
   return M
end

# Module-test.jl

function rand_homomorphism(M::AbstractAlgebra.FPModule{T}, vals...) where T <: RingElement
   rk = rand(1:5)
   m = ngens(M)
   R = base_ring(M)
   F = FreeModule(R, rk)
   S = MatrixSpace(R, rk, m)
   mat = rand(S, vals...)
   f = ModuleHomomorphism(F, M, mat)
   ngens1 = rand(1:3)
   gens1 = [rand(F, vals...) for j in 1:ngens1]
   S, g = sub(F, gens1)
   hom1 = compose(g, f)
   return S, hom1
end

# Map-test.jl

mutable struct MyMap <: Map{AbstractAlgebra.Integers{BigInt}, AbstractAlgebra.Integers{BigInt}, SetMap, MyMap}
   a::Int
end

Generic.domain(f::Map(MyMap)) = AbstractAlgebra.JuliaZZ
Generic.codomain(f::Map(MyMap)) = AbstractAlgebra.JuliaZZ

a(f::Map(MyMap)) = Generic.get_field(f, :a)

(f::MyMap)(x) =  a(f)*(x + 1)
