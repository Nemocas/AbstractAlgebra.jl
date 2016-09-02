###############################################################################
#
#   fmpz_mpoly.jl : Flint multivariate polynomials over fmpz
#
###############################################################################

export FmpzMPolyRing, fmpz_mpoly

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type{S, N}(::Type{fmpz_mpoly{S, N}}) = FmpzMPolyRing{S, N}

elem_type{S, N}(::FmpzMPolyRing{S, N}) = fmpz_mpoly{S, N}

vars(a::FmpzMPolyRing) = a.S

function gens{S, N}(R::FmpzMPolyRing{S, N})
   A = Array(fmpz_mpoly{S, N}, R.num_vars)
   for i = 1:R.num_vars
      z = R()
      ccall((:fmpz_mpoly_gen, :libflint), Void,
            (Ptr{fmpz_mpoly}, Int, Ptr{FmpzMPolyRing}), &z, i - 1, &R)
      A[i] = z
   end
   return A
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

num_vars(x::fmpz_mpoly) = parent(x).num_vars

###############################################################################
#
#   AbstractString{} I/O
#
###############################################################################

function show(io::IO, x::fmpz_mpoly)
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fmpz_mpoly_get_str_pretty, :libflint), Ptr{UInt8}, 
          (Ptr{fmpz_mpoly}, Ptr{Ptr{UInt8}}, Ptr{FmpzMPolyRing}), 
          &x, [bytestring(string(s)) for s in vars(parent(x))], &x.parent)

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{UInt8},), cstr)
   end
end

function show(io::IO, p::FmpzMPolyRing)
   const max_vars = 5 # largest number of variables to print
   n = p.num_vars
   print(io, "Multivariate Polynomial Ring in ")
   if n > max_vars
      print(io, p.num_vars)
      print(io, " variables ")
   end
   for i = 1:min(n - 1, max_vars - 1)
      print(io, string(p.S[i]), ", ")
   end
   if n > max_vars
      print(io, "..., ")
   end
   print(io, string(p.S[n]))
   print(io, " over ")
   show(io, base_ring(p))
end

###############################################################################
#
#   Basic arithmetic
#
###############################################################################

function +{S, N}(a::fmpz_mpoly{S, N}, b::fmpz_mpoly{S, N})
   z = parent(a)()
   ccall((:fmpz_mpoly_add, :libflint), Void, 
       (Ptr{fmpz_mpoly}, Ptr{fmpz_mpoly}, Ptr{fmpz_mpoly}, Ptr{FmpzMPolyRing}),
       &z, &a, &b, &a.parent)
   return z
end

function -{S, N}(a::fmpz_mpoly{S, N}, b::fmpz_mpoly{S, N})
   z = parent(a)()
   ccall((:fmpz_mpoly_sub, :libflint), Void, 
       (Ptr{fmpz_mpoly}, Ptr{fmpz_mpoly}, Ptr{fmpz_mpoly}, Ptr{FmpzMPolyRing}),
       &z, &a, &b, &a.parent)
   return z
end

###############################################################################
#
#   Ad hoc arithmetic
#
###############################################################################

function *(a::fmpz_mpoly, b::Int)
   r = parent(a)()
   ccall((:fmpz_mpoly_scalar_mul_si, :libflint), Void,
        (Ptr{fmpz_mpoly}, Ptr{fmpz_mpoly}, Int, Ptr{FmpzMPolyRing}), 
        &r, &a, b, &a.parent)
   return r
end

function *(a::fmpz_mpoly, b::fmpz)
   r = parent(a)()
   ccall((:fmpz_mpoly_scalar_mul_fmpz, :libflint), Void,
        (Ptr{fmpz_mpoly}, Ptr{fmpz_mpoly}, Ptr{fmpz}, Ptr{FmpzMPolyRing}), 
        &r, &a, &b, &a.parent)
   return r
end

*(a::fmpz, b::fmpz_mpoly) = b*a

*(a::Int, b::fmpz_mpoly) = b*a

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!{S, N}(a::fmpz_mpoly{S, N}, b::fmpz_mpoly{S, N})
   ccall((:fmpz_mpoly_add, :libflint), Void,
         (Ptr{fmpz_mpoly}, Ptr{fmpz_mpoly}, Ptr{fmpz_mpoly}),
                                                         &a, &a, &b, &a.parent)
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call{S, N}(R::FmpzMPolyRing{S, N})
   z = fmpz_mpoly{S, N}(R)
   z.parent = R
   return z
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

function PolynomialRing(R::FlintIntegerRing, s::Array{String, 1}; cached::Bool = true, ordering::Symbol = :lex)
   U = [Symbol(x) for x in s]
   N = (ordering == :deglex || ordering == :degrevlex) ? length(U) + 1 : length(U)
   # default to 8 bit exponent fields
   parent_obj = FmpzMPolyRing{ordering, N}(U, cached)
   return tuple(parent_obj, gens(parent_obj)...)
end



