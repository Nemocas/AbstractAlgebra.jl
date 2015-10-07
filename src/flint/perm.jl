###############################################################################
#
#   perm.jl : Flint permutation type
#
###############################################################################

export FlintPermGroup, perm, parity

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(a::perm) = a.parent

function deepcopy(a::perm)
   R = parent(a)
   p = R()
   ccall((:_perm_set, :libflint), Void, 
         (Ref{Int}, Ref{Int}, Int), p.d, a.d, R.n)
   return p
end

function parity(a::perm)
   R = parent(a)
   return Int(ccall((:_perm_parity, :libflint), Cint, 
                    (Ref{Int}, Int), a.d, R.n))
end

function getindex(a::perm, n::Int)
   return a.d[n] + 1
end
 
function setindex!(a::perm, d::Int, n::Int)
   a.d[n] = d - 1
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::FlintPermGroup)
   print(io, "Permutation group over ")
   print(io, R.n)
   print(io, " elements")
end

function show(io::IO, x::perm)
   print(io, "[")
   n = parent(x).n
   for i = 1:n
      print(io, x[i])
      if i != n
         print(io, ", ")
      end
   end
   print(io, "]")
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::perm, b::perm)
   R = parent(a)
   return Bool(ccall((:_perm_equal, :libflint), Cint, 
         (Ref{Int}, Ref{Int}, Int), a.d, b.d, R.n))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function *(a::perm, b::perm)
   R = parent(a)
   p = R()
   ccall((:_perm_compose, :libflint), Void, 
         (Ref{Int}, Ref{Int}, Ref{Int}, Int), p.d, a.d, b.d, R.n)
   return p
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::perm)
   R = parent(a)
   p = R()
   ccall((:_perm_inv, :libflint), Void, 
         (Ref{Int}, Ref{Int}, Int), p.d, a.d, R.n)
   return p
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function Base.call(R::FlintPermGroup)
   z = perm(R.n)
   z.parent = R
   return z
end

function Base.call(R::FlintPermGroup, a::Array{Int, 1})
   length(a) != R.n && error("Unable to coerce to permutation")
   z = perm(a)
   z.parent = R
   return z
end

function Base.call(R::FlintPermGroup, a::perm)
   parent(a) != R && error("Unable to coerce to permutation")
   return a
end

###############################################################################
#
#   FlintPermGroup constructor
#
###############################################################################

# handled by inner constructors
