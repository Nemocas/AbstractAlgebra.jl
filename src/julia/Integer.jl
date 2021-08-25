const JuliaZZ = Integers{BigInt}()

parent(a::T) where T <: Integer = Integers{T}()

elem_type(::Type{Integers{T}}) where T <: Integer = T

parent_type(::Type{T}) where T <: Integer = Integers{T}

base_ring(a::Integers{T}) where T <: Integer = Union{}

zero(::Integers{T}) where T <: Integer = T(0)

one(::Integers{T}) where T <: Integer = T(1)

canonical_unit(a::T) where T <: Integer = a < 0 ? T(-1) : T(1)

function divrem(a::BigInt, b::BigInt)
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q, r
end

function divrem(a::Int, b::Int)
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q, r
end

function divrem(a::BigInt, b::Int)
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q, r
end

function divrem(a::S, b::T) where {S <: Integer, T <: Integer}
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q, r
end

function div(a::S, b::T) where {S <: Integer, T <: Integer}
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q
end

function divexact(a::Integer, b::Integer; check::Bool=true)
   if check
      q, r = divrem(a, b)
      iszero(r) || throw(ArgumentError("Not an exact division"))
   else
      q = div(a, b)
   end
   return q
end

function divexact(a::BigInt, b::BigInt; check::Bool=true)
   q = BigInt()
   if check
      r = BigInt()
      ccall((:__gmpz_tdiv_qr, :libgmp), Nothing,
              (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), q, r, a, b)
      r != 0 && throw(ArgumentError("Not an exact division"))
   else
      ccall((:__gmpz_divexact, :libgmp), Nothing,
                              (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), q, a, b)
   end
   return q
end

function divexact(a::BigInt, b::Int; check::Bool=true)
   q = BigInt()
   sgn = b < 0
   if check
      r = BigInt()
      ccall((:__gmpz_tdiv_qr_ui, :libgmp), Nothing,
           (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Int), q, r, a, sgn ? -b : b)
      r != 0 && throw(ArgumentError("Not an exact division"))
   else
      ccall((:__gmpz_divexact_ui, :libgmp), Nothing,
                           (Ref{BigInt}, Ref{BigInt}, Int), q, a, sgn ? -b : b)
   end
   return sgn ? -q : q
end

function divexact(a::BigInt, b::UInt; check::Bool=true)
   q = BigInt()
   if check
      r = BigInt()
      ccall((:__gmpz_tdiv_qr_ui, :libgmp), Nothing,
                     (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, UInt), q, r, a, b)
      r != 0 && throw(ArgumentError("Not an exact division"))
   else
      ccall((:__gmpz_divexact_ui, :libgmp), Nothing,
                                     (Ref{BigInt}, Ref{BigInt}, UInt), q, a, b)
   end
   return q
end

function inv(a::T) where T <: Integer
   if a == 1
      return one(T)
   elseif a == -1
      return -one(T)
   end
   iszero(a) && throw(DivideError())
   throw(ArgumentError("not a unit"))
end

function gcdinv(a::T, b::T) where T <: Integer
   g, s, t = gcdx(a, b)
   return g, s
end

function zero!(a::Integer)
   return 0
end

function zero!(a::BigInt)
   ccall((:__gmpz_set_si, :libgmp), Nothing, (Ref{BigInt}, Int), a, 0)
   return a
end

function mul!(a::T, b::T, c::T) where T <: Integer
   return b*c
end

function mul!(a::BigInt, b::BigInt, c::BigInt)
   ccall((:__gmpz_mul, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

function add!(a::T, b::T, c::T) where T <: Integer
   return b + c
end

function add!(a::BigInt, b::BigInt, c::BigInt)
   ccall((:__gmpz_add, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

function addeq!(a::T, b::T) where T <: Integer
   return a + b
end

function addeq!(a::BigInt, b::BigInt)
   ccall((:__gmpz_add, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, a, b)
   return a
end

function addmul!(a::T, b::T, c::T, d::T) where T <: Integer
   return a + b*c
end

function addmul!(a::BigInt, b::BigInt, c::BigInt, d::BigInt)
   ccall((:__gmpz_addmul, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

function addmul!(a::T, b::T, c::T) where T <: Integer # special case, no temporary required
   return a + b*c
end

function addmul!(a::BigInt, b::BigInt, c::BigInt) # special case, no temporary required
   ccall((:__gmpz_addmul, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

function (a::Integers{T})() where T <: Integer
   return T(0)
end

function (a::Integers{T})(b::Integer) where T <: Integer
   return T(b)
end
