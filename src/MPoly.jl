###############################################################################
#
#   MPoly.jl : sparse distributed multivariate polynomials over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring(a::MPolyRingElem{T}) where T <: RingElement = base_ring(parent(a))

coefficient_ring(a::MPolyRingElem) = base_ring(a)

coefficient_ring(R::MPolyRing) = base_ring(R)

@doc raw"""
    mpoly_type(::Type{T}) where T<:RingElement
    mpoly_type(::T) where T<:RingElement
    mpoly_type(::Type{S}) where S<:Ring
    mpoly_type(::S) where S<:Ring

The type of multivariate polynomials with coefficients of type `T` respectively `elem_type(S)`.
Falls back to `Generic.MPoly{T}`.

See also [`mpoly_ring_type`](@ref), [`dense_poly_type`](@ref) and [`dense_poly_ring_type`](@ref).

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> mpoly_type(AbstractAlgebra.ZZ(1))
AbstractAlgebra.Generic.MPoly{BigInt}

julia> mpoly_type(elem_type(AbstractAlgebra.ZZ))
AbstractAlgebra.Generic.MPoly{BigInt}

julia> mpoly_type(AbstractAlgebra.ZZ)
AbstractAlgebra.Generic.MPoly{BigInt}

julia> mpoly_type(typeof(AbstractAlgebra.ZZ))
AbstractAlgebra.Generic.MPoly{BigInt}
```
"""
mpoly_type(::Type{T}) where T<:RingElement = Generic.MPoly{T}
mpoly_type(::Type{S}) where S<:Ring = mpoly_type(elem_type(S))
mpoly_type(x) = mpoly_type(typeof(x)) # to stop this method from eternally recursing on itself, we better add ...
mpoly_type(::Type{T}) where T = throw(ArgumentError("Type `$T` must be subtype of `RingElement`."))

@doc raw"""
    mpoly_ring_type(::Type{T}) where T<:RingElement
    mpoly_ring_type(::T) where T<:RingElement
    mpoly_ring_type(::Type{S}) where S<:Ring
    mpoly_ring_type(::S) where S<:Ring

The type of multivariate polynomial rings with coefficients of type `T`
respectively `elem_type(S)`. Implemented via [`mpoly_type`](@ref).

See also [`dense_poly_type`](@ref) and [`dense_poly_ring_type`](@ref).

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> mpoly_ring_type(AbstractAlgebra.ZZ(1))
AbstractAlgebra.Generic.MPolyRing{BigInt}

julia> mpoly_ring_type(elem_type(AbstractAlgebra.ZZ))
AbstractAlgebra.Generic.MPolyRing{BigInt}

julia> mpoly_ring_type(AbstractAlgebra.ZZ)
AbstractAlgebra.Generic.MPolyRing{BigInt}

julia> mpoly_ring_type(typeof(AbstractAlgebra.ZZ))
AbstractAlgebra.Generic.MPolyRing{BigInt}
```
"""
mpoly_ring_type(x) = parent_type(mpoly_type(x))

function is_domain_type(::Type{T}) where {S <: RingElement, T <: AbstractAlgebra.MPolyRingElem{S}}
   return is_domain_type(S)
end

function is_exact_type(a::Type{T}) where {S <: RingElement, T <: AbstractAlgebra.MPolyRingElem{S}}
   return is_exact_type(S)
end

@doc raw"""
    vars(p::AbstractAlgebra.MPolyRingElem{T}) where {T <: RingElement}

Return the variables actually occurring in $p$.
"""
function vars(p::AbstractAlgebra.MPolyRingElem{T}) where {T <: RingElement}
   U = typeof(p)
   vars_in_p = Vector{U}(undef, 0)
   n = nvars(p.parent)
   gen_list = gens(p.parent)
   biggest = [0 for i in 1:n]
   for v in exponent_vectors(p)
      for j = 1:n
         if v[j] > biggest[j]
            biggest[j] = v[j]
         end
      end
   end
   for i = 1:n
      if biggest[i] != 0
         push!(vars_in_p, gen_list[i])
      end
   end
   return(vars_in_p)
end

@doc raw"""
    var_index(p::AbstractAlgebra.MPolyRingElem{T}) where {T <: RingElement}

Return the index of the given variable $x$. If $x$ is not a variable
in a multivariate polynomial ring, an exception is raised.
"""
function var_index(x::AbstractAlgebra.MPolyRingElem{T}) where {T <: RingElement}
   !is_monomial(x) && error("Not a variable in var_index")
   exps = first(exponent_vectors(x))
   count = 0
   index = 0
   for i = 1:length(exps)
      if exps[i] > 1
         error("Not a variable in var_index")
      end
      if exps[i] == 1
         count += 1
         index = i
      end
   end
   count != 1 && error("Not a variable in var_index")
   return index
end

function characteristic(a::MPolyRing{T}) where T <: RingElement
   return characteristic(base_ring(a))
end

###############################################################################
#
#   Multivariate coefficients
#
###############################################################################

@doc raw"""
    coeff(a::AbstractAlgebra.MPolyRingElem{T}, vars::Vector{Int}, exps::Vector{Int}) where T <: RingElement

Return the "coefficient" of $a$ (as a multivariate polynomial in the same
ring) of the monomial consisting of the product of the variables of the given
indices raised to the given exponents (note that not all variables need to
appear and the exponents can be zero). E.g. `coeff(f, [1, 3], [0, 2])` returns
the coefficient of $x^0*z^2$ in the polynomial $f$ (assuming variables
$x, y, z$ in that order).
"""
function coeff(a::AbstractAlgebra.MPolyRingElem{T}, vars::Vector{Int}, exps::Vector{Int}) where T <: RingElement
   unique(vars) != vars && error("Variables not unique")
   length(vars) != length(exps) &&
       error("Number of variables does not match number of exponents")
   for i = 1:length(vars)
      if vars[i] < 1 || vars[i] > nvars(parent(a))
         error("Variable index not in range")
      end
      if exps[i] < 0
         error("Exponent cannot be negative")
      end
   end
   S = parent(a)
   M = Generic.MPolyBuildCtx(S)
   cvzip = zip(coefficients(a), exponent_vectors(a))
   for (c, v) in cvzip
      flag = true
      for j = 1:length(vars)
         if v[vars[j]] != exps[j]
            flag = false
            break
         else
            v[vars[j]] = 0
         end
      end
      if flag
         push_term!(M, c, v)
      end
   end
   return finish(M)
end

@doc raw"""
    coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where T <: AbstractAlgebra.MPolyRingElem

Return the "coefficient" of $a$ (as a multivariate polynomial in the same
ring) of the monomial consisting of the product of the given variables
to the given exponents (note that not all variables need to appear and the
exponents can be zero). E.g. `coeff(f, [x, z], [0, 2])` returns the
coefficient of $x^0*z^2$ in the polynomial $f$.
"""
function coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where T <: AbstractAlgebra.MPolyRingElem
   varidx = [var_index(x) for x in vars]
   return coeff(a, varidx, exps)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

# Fallback hash function for multivariate polynomials implementing the
# iterators
function Base.hash(x::MPolyRingElem{T}, h::UInt) where {T <: RingElement}
   b = 0x53dd43cd511044d1%UInt
   for (e, c) in zip(exponent_vectors(x), coefficients(x))
      b = xor(b, hash(c, h), h)
      b = xor(b, hash(e, h), h)
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

@doc raw"""
    is_degree(s::Symbol)

Return `true` if the given symbol represents a degree ordering (deglex or
degrevlex).
"""
is_degree(s::Symbol) = s == :deglex || s == :degrevlex

@doc raw"""
    is_reverse(s::Symbol)

Return `true` if the given symbol represents a reverse ordering (degrevlex).
"""
is_reverse(s::Symbol) = s == :degrevlex

@doc raw"""
    coeff(f::AbstractAlgebra.MPolyRingElem{T}, m::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Return the coefficient of the monomial $m$ of the polynomial $f$. If there
is no such monomial, zero is returned.
"""
function coeff(f::AbstractAlgebra.MPolyRingElem{T}, m::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
    !is_monomial(m) && error("Not a monomial in coeff")
    v1 = first(exponent_vectors(m))
    cvzip = zip(coefficients(f), exponent_vectors(f))
    for (c, v) in cvzip
        if v == v1
            return c
        end
    end
    return zero(base_ring(f))
end

@doc raw"""
    leading_coefficient(p::MPolyRingElem)

Return the leading coefficient of the polynomial $p$.
"""
function leading_coefficient(p::MPolyRingElem{T}) where T <: RingElement
   if iszero(p)
      return zero(base_ring(p))
   else
      return first(coefficients(p))
   end
end

@doc raw"""
    trailing_coefficient(p::MPolyRingElem)

Return the trailing coefficient of the polynomial $p$, i.e. the coefficient of
the last nonzero term, or zero if the polynomial is zero.
"""
function trailing_coefficient(p::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   coeff = zero(base_ring(p))
   for c in coefficients(p)
      coeff = c
   end
   return coeff
end

@doc raw"""
    tail(p::MPolyRingElem)

Return the tail of the polynomial $p$, i.e. the polynomial without its leading
term (if any).
"""
function tail(p::MPolyRingElem{T}) where T <: RingElement
   S = parent(p)
   if iszero(p)
      return S()
   end
   ctx = Generic.MPolyBuildCtx(S)
   tail_cv = Iterators.drop(zip(coefficients(p), exponent_vectors(p)), 1)
   for (c, v) in tail_cv
      push_term!(ctx, c, v)
   end
   return finish(ctx)
end

@doc raw"""
    constant_coefficient(p::MPolyRingElem)

Return the constant coefficient of the polynomial $p$ or zero if it doesn't
have one.
"""
function constant_coefficient(p::MPolyRingElem{T}) where T <: RingElement
   if !iszero(p)
      for (c, v) in zip(coefficients(p), exponent_vectors(p))
         if iszero(v)
            return c
         end
      end
   end
   return zero(base_ring(p))
end

function constant_coefficient(p::MPolyRingElem)
   len = length(p)
   if !iszero(p) && iszero(exponent_vector(p, len))
      return coeff(p, len)
   end
   return zero(base_ring(p))
end

@doc raw"""
    leading_monomial(p::MPolyRingElem)

Return the leading monomial of $p$.
This function throws an `ArgumentError` if $p$ is zero.
"""
function leading_monomial(p::MPolyRingElem{T}) where T <: RingElement
   if iszero(p)
      throw(ArgumentError("Zero polynomial does not have a leading monomial"))
   end
   return first(monomials(p))
end

@doc raw"""
    leading_exponent_vector(p::MPolyRingElem)

Return the exponent vector of the leading term of $p$. The return is a Julia
1-dimensional array giving the exponent for each variable of the leading term.
This function throws an `ArgumentError` if $p$ is zero.
"""
function leading_exponent_vector(p::MPolyRingElem{T}) where T <: RingElement
   if iszero(p)
      throw(ArgumentError("Zero polynomial does not have a leading exponent vector"))
   end
   return first(exponent_vectors(p))
end

@doc raw"""
    leading_term(p::MPolyRingElem)

Return the leading term of the polynomial p.
This function throws an `ArgumentError` if $p$ is zero.
"""
function leading_term(p::MPolyRingElem{T}) where T <: RingElement
   if iszero(p)
      throw(ArgumentError("Zero polynomial does not have a leading term"))
   end
   return first(terms(p))
end

@doc raw"""
    degree(f::AbstractAlgebra.MPolyRingElem{T}, i::Int) where T <: RingElement

Return the degree of the polynomial $f$ in terms of the i-th variable.
"""
function degree(f::AbstractAlgebra.MPolyRingElem{T}, i::Int) where T <: RingElement
   biggest = -1
   if length(f) != 0
      R = parent(f)
      if ordering(R) == :lex && i == 1
         biggest = first(exponent_vectors(f))[1]
      else
         for v in exponent_vectors(f)
            if v[i] > biggest
               biggest = v[i]
            end
         end
      end
   end
   return biggest
end

@doc raw"""
    degree(f::AbstractAlgebra.MPolyRingElem{T}, x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Return the degree of the polynomial $f$ in terms of the variable $x$.
"""
function degree(f::AbstractAlgebra.MPolyRingElem{T}, x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   return degree(f, var_index(x))
end

@doc raw"""
    degrees(f::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Return an array of the degrees of the polynomial $f$ in terms of each variable.
"""
function degrees(f::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   R = parent(f)
   if nvars(R) == 1 && ordering(R) == :lex && length(f) > 0
      return first(exponent_vectors(f))
   else
      biggest = [-1 for i = 1:nvars(R)]
      for v in exponent_vectors(f)
         for j = 1:nvars(R)
            if v[j] > biggest[j]
               biggest[j] = v[j]
            end
         end
      end
      return biggest
   end
end

one(R::MPolyRing) = R(1)

zero(R::MPolyRing) = R(0)

function isone(x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   return length(x) == 1 && iszero(first(exponent_vectors(x))) &&
          first(coefficients(x)) == 1
end

iszero(x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement = length(x) == 0

function is_unit(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   if is_constant(a)
      return is_unit(leading_coefficient(a))
   elseif is_domain_type(elem_type(coefficient_ring(a)))
      return false
   elseif length(a) == 1
      return false
   else
      throw(NotImplementedError(:is_unit, a))
   end
end

function content(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   z = zero(coefficient_ring(a))
   for c in coefficients(a)
      z = gcd(z, c)
   end
   return z
end

function is_zero_divisor(x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   return is_zero_divisor(content(x))
end

function is_zero_divisor_with_annihilator(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   f, b = is_zero_divisor_with_annihilator(content(a))
   return f, parent(a)(b)
end

@doc raw"""
    is_constant(x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Return `true` if `x` is a degree zero polynomial or the zero polynomial, i.e.
a constant polynomial.
"""
function is_constant(x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
    return length(x) == 0 || (length(x) == 1 &&
                        iszero(first(exponent_vectors(x))))
end

@doc raw"""
    is_term(x::MPoly)

Return `true` if the given polynomial has precisely one term.
"""
is_term(x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement = length(x) == 1

@doc raw"""
    is_monomial(x::AbstractAlgebra.MPolyRingElem)

Return `true` if the given polynomial has precisely one term whose coefficient is one.
"""
function is_monomial(x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   return length(x) == 1 && isone(first(coefficients(x)))
end

###############################################################################
#
#   Iterators
#
###############################################################################

@doc raw"""
    coefficients(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Return an iterator for the coefficients of the given polynomial. To retrieve
an array of the coefficients, use `collect(coefficients(a))`.
"""
function coefficients(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   return Generic.MPolyCoeffs(a)
end

@doc raw"""
    exponent_vectors(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Return an iterator for the exponent vectors of the given polynomial. To
retrieve an array of the exponent vectors, use
`collect(exponent_vectors(a))`.
"""
function exponent_vectors(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   return Generic.MPolyExponentVectors(a)
end

@doc raw"""
    monomials(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Return an iterator for the monomials of the given polynomial. To retrieve
an array of the monomials, use `collect(monomials(a))`.
"""
function monomials(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   return Generic.MPolyMonomials(a)
end

@doc raw"""
    terms(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Return an iterator for the terms of the given polynomial. To retrieve
an array of the terms, use `collect(terms(a))`.
"""
function terms(a::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   return Generic.MPolyTerms(a)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::MPolyRingElem, x = symbols(parent(a)); context = nothing)
   sum = Expr(:call, :+)
   n = nvars(parent(a))
   for (c, v) in zip(coefficients(a), exponent_vectors(a))
      prod = Expr(:call, :*)
      if !isone(c)
         push!(prod.args, expressify(c, context = context))
      end
      for i in 1:n
         if v[i] > 1
            push!(prod.args, Expr(:call, :^, x[i], v[i]))
         elseif v[i] == 1
            push!(prod.args, x[i])
         end
      end
      push!(sum.args, prod)
   end
   return sum
end

@enable_all_show_via_expressify MPolyRingElem


function show(io::IO, p::MPolyRing)
   local max_vars = 5 # largest number of variables to print
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
   print(IOContext(io, :compact => true), base_ring(p))
end

function canonical_unit(x::AbstractAlgebra.MPolyRingElem)
   if length(x) == 0
      return base_ring(x)()
   else
      return canonical_unit(coeff(x, 1))
   end
end

###############################################################################
#
#   Inflation/deflation
#
###############################################################################

@doc raw"""
    deflation(f::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Compute deflation parameters for the exponents of the polynomial $f$. This
is a pair of arrays of integers, the first array of which (the shift) gives
the minimum exponent for each variable of the polynomial, and the second of
which (the deflation) gives the gcds of all the exponents after subtracting
the shift, again per variable. This functionality is used by gcd (and can be
used by factorisation algorithms).
"""
function deflation(f::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   N = nvars(parent(f))
   if length(f) == 0
      return [0 for i in 1:N], [0 for i in 1:N]
   end
   defl = [0 for i in 1:N]
   shift = first(exponent_vectors(f))
   for v in Iterators.drop(exponent_vectors(f), 1)
      for j = 1:N
         exj = v[j]
         if exj < shift[j]
            defl[j] = defl[j] == 1 ? 1 : gcd(defl[j], shift[j] - exj)
            shift[j] = exj
         else
            defl[j] = defl[j] == 1 ? 1 : gcd(defl[j], exj - shift[j])
         end
      end
   end
   return shift, defl
end

@doc raw"""
    deflate(f::AbstractAlgebra.MPolyRingElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement

Return a polynomial with the same coefficients as $f$ but whose exponents
have been reduced by the given shifts (supplied as an array of shifts,
one for each variable), then deflated (divided) by the given exponents
(again supplied as an array of deflation factors, one for each variable).
The algorithm automatically replaces a deflation of $0$ by $1$, to avoid
division by $0$.
"""
function deflate(f::AbstractAlgebra.MPolyRingElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement
   S = parent(f)
   N = nvars(S)
   for i = 1:N
      if defl[i] == 0
         defl[i] = 1
      end
   end
   M = Generic.MPolyBuildCtx(S)
   cvzip = zip(coefficients(f), exponent_vectors(f))
   for (c, v) in cvzip
      for j = 1:N
         v[j] = div(v[j] - shift[j], defl[j])
      end
      push_term!(M, c, v)
   end
   return finish(M)
end

@doc raw"""
    deflate(f::AbstractAlgebra.MPolyRingElem{T}, defl::Vector{Int}) where T <: RingElement

Return a polynomial with the same coefficients as $f$ but whose exponents have
been deflated (divided) by the given exponents (supplied as an array of
deflation factors, one for each variable).

The algorithm automatically replaces a deflation of $0$ by $1$, to avoid
division by $0$.
"""
function deflate(f::AbstractAlgebra.MPolyRingElem{T}, defl::Vector{Int}) where T <: RingElement
   return deflate(f, [0 for i in 1:nvars(parent(f))], defl)
end

@doc raw"""
    deflate(f::AbstractAlgebra.MPolyRingElem{T}, defl::Vector{Int}) where T <: RingElement

Return a polynomial with the same coefficients as $f$ but whose exponents have
been deflated maximally, i.e. with each exponent divide by the largest integer
which divides the degrees of all exponents of that variable in $f$.
"""
function deflate(f::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   shift, defl = deflation(f)
   defl = gcd.(shift, defl)
   return deflate(f, defl), defl
end

function inflate_deflate_vectors(R::AbstractAlgebra.MPolyRing, vars::Vector{Int}, shift::Vector{Int}, defl::Vector{Int})::Tuple{Vector{Int},Vector{Int}}
   unique(vars) != vars && error("Variables not unique")
   !(length(vars) == length(shift) == length(defl)) && error("Number of variables does not match lengths of shift and deflation vectors")
   shift1 = zeros(Int, nvars(R))
   defl1 = ones(Int, nvars(R))
   for i in 1:(length(vars))
      !(1 <= vars[i] <= nvars(R)) && error("Variable index not in range")
      shift1[vars[i]] = shift[i]
      defl1[vars[i]] = defl[i]
   end
   return (shift1, defl1)
end

@doc raw"""
    deflate(f::AbstractAlgebra.MPolyRingElem, vars::Vector{Int}, shift::Vector{Int}, defl::Vector{Int})

Return a polynomial with the same coefficients as $f$ but where exponents of
some variables (supplied as an array of variable indices) have been reduced by
the given shifts (supplied as an array of shifts), then deflated (divided) by
the given exponents (again supplied as an array of deflation factors). The
algorithm automatically replaces a deflation of $0$ by $1$, to avoid division by
$0$.
"""
function deflate(f::AbstractAlgebra.MPolyRingElem, vars::Vector{Int}, shift::Vector{Int}, defl::Vector{Int})
   (shift1, defl1) = inflate_deflate_vectors(parent(f), vars, shift, defl)
   return deflate(f, shift1, defl1)
end

@doc raw"""
    deflate(f::T, vars::Vector{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: AbstractAlgebra.MPolyRingElem

Return a polynomial with the same coefficients as $f$ but where the exponents of
the given variables have been reduced by the given shifts (supplied as an array
of shifts), then deflated (divided) by the given exponents (again supplied as an
array of deflation factors). The algorithm automatically replaces a deflation of
$0$ by $1$, to avoid division by $0$.
"""
function deflate(f::T, vars::Vector{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: AbstractAlgebra.MPolyRingElem
   varidx = [var_index(x) for x in vars]
   return deflate(f, varidx, shift, defl)
end


@doc raw"""
    inflate(f::AbstractAlgebra.MPolyRingElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement

Return a polynomial with the same coefficients as $f$ but whose exponents
have been inflated (multiplied) by the given deflation exponents (supplied
as an array of inflation factors, one for each variable) and then increased
by the given shifts (again supplied as an array of shifts, one for each
variable).
"""
function inflate(f::AbstractAlgebra.MPolyRingElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement
   S = parent(f)
   N = nvars(S)
   M = Generic.MPolyBuildCtx(S)
   cvzip = zip(coefficients(f), exponent_vectors(f))
   for (c, v) in cvzip
      for j = 1:N
         v[j] = v[j]*defl[j] + shift[j]
      end
      push_term!(M, c, v)
   end
   return finish(M)
end

@doc raw"""
    inflate(f::AbstractAlgebra.MPolyRingElem{T}, defl::Vector{Int}) where T <: RingElement

Return a polynomial with the same coefficients as $f$ but whose exponents
have been inflated (multiplied) by the given deflation exponents (supplied
as an array of inflation factors, one for each variable).
"""
function inflate(f::AbstractAlgebra.MPolyRingElem{T}, defl::Vector{Int}) where T <: RingElement
   return inflate(f, [0 for i in 1:nvars(parent(f))], defl)
end

@doc raw"""
    inflate(f::AbstractAlgebra.MPolyRingElem, vars::Vector{Int}, shift::Vector{Int}, defl::Vector{Int})

Return a polynomial with the same coefficients as $f$ but where exponents of
some variables (supplied as an array of variable indices) have been inflated
(multiplied) by the given deflation exponents (supplied as an array of inflation
factors) and then increased by the given shifts (again supplied as an array of
shifts).
"""
function inflate(f::AbstractAlgebra.MPolyRingElem, vars::Vector{Int}, shift::Vector{Int}, defl::Vector{Int})
   (shift1, defl1) = inflate_deflate_vectors(parent(f), vars, shift, defl)
   return inflate(f, shift1, defl1)
end

@doc raw"""
    inflate(f::T, vars::Vector{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: AbstractAlgebra.MPolyRingElem

Return a polynomial with the same coefficients as $f$ but where the exponents of
the given variables have been inflated (multiplied) by the given deflation
exponents (supplied as an array of inflation factors) and then increased by the
given shifts (again supplied as an array of shifts).
"""
function inflate(f::T, vars::Vector{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: AbstractAlgebra.MPolyRingElem
   varidx = [var_index(x) for x in vars]
   return inflate(f, varidx, shift, defl)
end


################################################################################
#
#   Remove and valuation
#
################################################################################

@doc raw"""
    remove(z::MPolyRingElem{T}, p::MPolyRingElem{T}) where {T <: RingElement}

Compute the valuation of $z$ at $p$, that is, the largest $k$ such that
$p^k$ divides $z$. Additionally, $z/p^k$ is returned as the second return.
value.

See also `valuation`, which only returns the valuation.
"""
function remove(z::MPolyRingElem{T}, p::MPolyRingElem{T}) where {T <: RingElement}
   check_parent(z, p)
   iszero(z) && error("Not yet implemented")
   fl, q = divides(z, p)
   if !fl
      return 0, z
   end
   v = 0
   qn = q
   while fl
      q = qn
      fl, qn = divides(q, p)
      v += 1
   end
   return v, q
end

@doc raw"""
    valuation(z::MPolyRingElem{T}, p::MPolyRingElem{T}) where {T <: RingElement}

Compute the valuation of $z$ at $p$, that is, the largest $k$ such that
$p^k$ divides $z$.

See also `remove`, which also returns $z/p^k$.
"""
function valuation(z::MPolyRingElem{T}, p::MPolyRingElem{T}) where {T <: RingElement}
  v, _ = remove(z, p)
  return v
end

###############################################################################
#
#   Evaluation
#
###############################################################################

@doc raw"""
    evaluate(a::AbstractAlgebra.MPolyRingElem{T}, vals::Vector{U}) where {T <: RingElement, U <: RingElement}

Evaluate the polynomial expression by substituting in the array of values for
each of the variables. The evaluation will succeed if multiplication is
defined between elements of the coefficient ring of $a$ and elements of the
supplied vector.
"""
function evaluate(a::AbstractAlgebra.MPolyRingElem{T}, vals::Vector{U}) where {T <: RingElement, U <: RingElement}
   length(vals) != nvars(parent(a)) && error("Incorrect number of values in evaluation")
   R = base_ring(a)
   if (U <: Integer && U != BigInt) ||
      (U <: Rational && U != Rational{BigInt})
      c = zero(R)*zero(U)
      V = typeof(c)
      if U != V
         vals = [parent(c)(v) for v in vals]
         powers = [Dict{Int, V}() for i in 1:length(vals)]
      else
         powers = [Dict{Int, U}() for i in 1:length(vals)]
      end
   else
      powers = [Dict{Int, U}() for i in 1:length(vals)]
   end
   # The best we can do here is to cache previously used powers of the values
   # being substituted, as we cannot assume anything about the relative
   # performance of powering vs multiplication. The function should not try
   # to optimise computing new powers in any way.
   # Note that this function accepts values in a non-commutative ring, so operations
   # must be done in a certain order.
   # But addition is associative.
   S = parent(one(R)*one(parent(vals[1])))
   r = elem_type(S)[zero(S)]
   i = UInt(1)
   cvzip = zip(coefficients(a), exponent_vectors(a))
   for (c, v) in cvzip
      t = one(S)
      for j = 1:length(vals)
         exp = v[j]
         if iszero(exp)
           continue
         end
         if !haskey(powers[j], exp)
            powers[j][exp] = vals[j]^exp
         end
         t = t*powers[j][exp]
      end
      push!(r, c*t)
      j = i = i + 1
      while iseven(j) && length(r) > 1
          top = pop!(r)
          r[end] = addeq!(r[end], top)
          j >>= 1
      end
   end
   while length(r) > 1
      top = pop!(r)
      r[end] = addeq!(r[end], top)
   end
   return r[1]
end

@doc raw"""
    evaluate(a::AbstractAlgebra.MPolyRingElem{T}, vars::Vector{Int}, vals::Vector{U}) where {T <: RingElement, U <: RingElement}

Evaluate the polynomial expression by substituting in the supplied values in
the array `vals` for the corresponding variables with indices given by the
array `vars`. The evaluation will succeed if multiplication is defined between
elements of the coefficient ring of $a$ and elements of `vals`.
"""
function evaluate(a::AbstractAlgebra.MPolyRingElem{T}, vars::Vector{Int}, vals::Vector{U}) where {T <: RingElement, U <: RingElement}
   unique(vars) != vars && error("Variables not unique")
   length(vars) != length(vals) &&
      error("Number of variables does not match number of values")
   for i = 1:length(vars)
      if vars[i] < 1 || vars[i] > nvars(parent(a))
         error("Variable index not in range")
      end
   end

   if length(vars) == 0
     return a
   end

   S = parent(a)
   R = base_ring(a)
   return _evaluate(a, S, R, vars, vals)
end

function _evaluate(a, S, R, vars, vals::Vector{U}) where {U <: Integer}
  c = zero(R) * zero(U)
  V = typeof(c)
  if V === U
     powers = Dict{Int, U}[Dict{Int, U}() for i in 1:length(vals)]
     return __evaluate(a, vars, vals, powers)
  else
     vals2 = V[parent(c)(v) for v in vals]
     powers = Dict{Int, V}[Dict{Int, V}() for i in 1:length(vals)]
     return __evaluate(a, vars, vals2, powers)
  end
end

function _evaluate(a, S, R, vars, vals::Vector{U}) where {U <: Rational}
  c = zero(R) * zero(U)
  V = typeof(c)
  if V === U
     powers = Dict{Int, U}[Dict{Int, U}() for i in 1:length(vals)]
     return __evaluate(a, vars, vals, powers)
  else
     vals2 = V[parent(c)(v) for v in vals]
     powers = Dict{Int, V}[Dict{Int, V}() for i in 1:length(vals)]
     return __evaluate(a, vars, vals2, powers)
  end
end

function _evaluate(a, S, R, vars, vals::Vector{U}) where {U <: RingElement}
  powers = Dict{Int, U}[Dict{Int, U}() for i in 1:length(vals)]
  return __evaluate(a, vars, vals, powers)
end

function _evaluate(a, S, R, vars, vals::Vector{BigInt})
  powers = Dict{Int, BigInt}[Dict{Int, BigInt}() for i in 1:length(vals)]
  return __evaluate(a, vars, vals, powers)
end

function _evaluate(a, S, R, vars, vals::Vector{Rational{BigInt}})
  powers = Dict{Int, Rational{BigInt}}[Dict{Int, Rational{BigInt}}() for i in 1:length(vals)]
  return __evaluate(a, vars, vals, powers)
end

function __evaluate(a, vars, vals, powers)
   R = base_ring(a)
   S = parent(a)
   # The best we can do here is to cache previously used powers of the values
   # being substituted, as we cannot assume anything about the relative
   # performance of powering vs multiplication. The function should not try
   # to optimise computing new powers in any way.

   # We use a geobucket if the result will be an element in the same ring as a
   if parent(vals[1] * one(S)) == S
     r = Generic.geobucket(S)
     cvzip = zip(coefficients(a), exponent_vectors(a))
     for (c, v) in cvzip
        t = one(S)
        for j = 1:length(vars)
           varnum = vars[j]
           exp = v[varnum]
           if iszero(exp)
             v[varnum] = 0
             continue
           end
           if !haskey(powers[j], exp)
              powers[j][exp] = vals[j]^exp
           end
           t *= powers[j][exp]
           v[varnum] = 0
        end
        M = Generic.MPolyBuildCtx(S)
        push_term!(M, c, v)
        push!(r, t*finish(M))
     end
     return finish(r)
   else
     K = parent(one(S) * vals[1])
     r = zero(K)
     cvzip = zip(coefficients(a), exponent_vectors(a))
     for (c, v) in cvzip
        t = one(K)
        for j = 1:length(vars)
           varnum = vars[j]
           exp = v[varnum]
           if iszero(exp)
             v[varnum] = 0
             continue
           end
           if !haskey(powers[j], exp)
              powers[j][exp] = vals[j]^exp
           end
           t *= powers[j][exp]
           v[varnum] = 0
        end
        M = Generic.MPolyBuildCtx(S)
        push_term!(M, c, v)
        addeq!(r, t*finish(M))
     end
     return r
   end
end

@doc raw"""
    evaluate(a::S, vars::Vector{S}, vals::Vector{U}) where {S <: AbstractAlgebra.MPolyRingElem{T}, U <: RingElement} where T <: RingElement

Evaluate the polynomial expression by substituting in the supplied values in
the array `vals` for the corresponding variables (supplied as polynomials)
given by the array `vars`. The evaluation will succeed if multiplication is
defined between elements of the coefficient ring of $a$ and elements of `vals`.
"""
function evaluate(a::S, vars::Vector{S}, vals::Vector{U}) where {S <: AbstractAlgebra.MPolyRingElem{T}, U <: RingElement} where T <: RingElement
   varidx = Int[var_index(x) for x in vars]
   return evaluate(a, varidx, vals)
end

@doc raw"""
    evaluate(a::AbstractAlgebra.MPolyRingElem{T}, vals::Vector{U}) where {T <: RingElement, U <: NCRingElem}

Evaluate the polynomial expression at the supplied values, which may be any
ring elements, commutative or non-commutative, but in the same ring. Evaluation
always proceeds in the order of the variables as supplied when creating the
polynomial ring to which $a$ belongs. The evaluation will succeed if a product
of a coefficient of the polynomial by one of the values is defined.
"""
function evaluate(a::AbstractAlgebra.MPolyRingElem{T}, vals::Vector{U}) where {T <: RingElement, U <: NCRingElem}
   return a(vals...)
end

################################################################################
#
#  Derivative
#
################################################################################

@doc raw"""
    derivative(f::AbstractAlgebra.MPolyRingElem{T}, j::Int) where {T <: RingElement}

Return the partial derivative of `f` with respect to $j$-th variable
of the polynomial ring.
"""
function derivative(f::AbstractAlgebra.MPolyRingElem{T}, j::Int) where T <: RingElement
   R = parent(f)
   iterz = zip(coefficients(f), exponent_vectors(f))
   Ctx = Generic.MPolyBuildCtx(R)
   for (c, v) in iterz
      if v[j] >= 1
         prod = c*v[j]
         v[j] -= 1
         push_term!(Ctx, prod, v)
      end
   end
   return finish(Ctx)
end

@doc raw"""
    derivative(f::AbstractAlgebra.MPolyRingElem{T}, x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement

Return the partial derivative of `f` with respect to `x`. The value `x` must
be a generator of the polynomial ring of `f`.
"""
function derivative(f::AbstractAlgebra.MPolyRingElem{T}, x::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   return derivative(f, var_index(x))
end


###############################################################################
#
#   Univariate polynomials
#
###############################################################################

@doc raw"""
    to_univariate(R::AbstractAlgebra.PolyRing{T}, p::AbstractAlgebra.MPolyRingElem{T}) where T <: AbstractAlgebra.RingElement

Assuming the polynomial $p$ is actually a univariate polynomial, convert the
polynomial to a univariate polynomial in the given univariate polynomial ring
$R$. An exception is raised if the polynomial $p$ involves more than one
variable.
"""
function to_univariate(R::AbstractAlgebra.PolyRing{T}, p::AbstractAlgebra.MPolyRingElem{T}) where T <: AbstractAlgebra.RingElement
   if !is_univariate(p)
      error("Can only convert univariate polynomials of type MPoly.")
   end
   if is_constant(p)
      return R(leading_coefficient(p))
   end
   return R(coefficients_of_univariate(p))
end

@doc raw"""
    is_univariate(p::AbstractAlgebra.MPolyRingElem)

Returns `true` if $p$ is a univariate polynomial, i.e. involves at most one
variable (thus constant polynomials are considered univariate), and `false`
otherwise. The result depends on the terms of the polynomial, not simply on
the number of variables in the polynomial ring.
"""
function is_univariate(p::AbstractAlgebra.MPolyRingElem{T}) where T <: RingElement
   if is_constant(p)
      return true
   end
   var = -1
   for v in exponent_vectors(p)
      n = count(x -> x != 0, v)
      if n > 1
         return false
      elseif n == 1
         if var == -1
            var = findfirst(x -> x != 0, v)
         elseif v[var] == 0
            return false
         end
      end
   end
   return true
end

@doc raw"""
    is_univariate(R::AbstractAlgebra.MPolyRing)

Returns `true` if $R$ is a univariate polynomial ring, i.e. has exactly one
variable, and `false` otherwise.
"""
function is_univariate(R::AbstractAlgebra.MPolyRing{T}) where T <: RingElement
   return nvars(R) == 1
end

@doc raw"""
    coefficients_of_univariate(p::MPolyRingElem)

Return the coefficients of p, which is assumed to be univariate,
as an array in ascending order.
"""
function coefficients_of_univariate(p::AbstractAlgebra.MPolyRingElem,
                                    check_univariate::Bool=true)

   coeffs = Vector{elem_type(coefficient_ring(p))}(undef, 0)
   var_index = -1
   for (c, v) in zip(coefficients(p), exponent_vectors(p))
      e = 0
      if var_index < 0
         for i in 1:length(v)
            if v[i] != 0
               e = v[i]
               var_index = i
               break
            end
         end
      else
         e = v[var_index]
      end

      while length(coeffs) <= e
         push!(coeffs, zero(coefficient_ring(p)))
      end
      coeffs[1 + e] = c

      if check_univariate
         for i in 1:length(v)
            if i != var_index && v[i] != 0
               error("Polynomial is not univariate.")
            end
         end
      end
   end
   return coeffs
end

################################################################################
#
#  Change base ring
#
################################################################################

function _change_mpoly_ring(R, Rx, cached)
   P, _ = AbstractAlgebra.polynomial_ring(R, map(string, symbols(Rx)), ordering = ordering(Rx), cached = cached)
   return P
end

@doc raw"""
    change_base_ring(R::Ring, p::MPolyRingElem{<: RingElement}; parent::MPolyRing, cached::Bool=true)

Return the polynomial obtained by coercing the non-zero coefficients of `p`
into `R`.

If the optional `parent` keyword is provided, the polynomial will be an
element of `parent`. The caching of the parent object can be controlled
via the `cached` keyword argument.
"""
function change_base_ring(R::Ring, p::MPolyRingElem{T}; cached::Bool=true, parent::AbstractAlgebra.MPolyRing = _change_mpoly_ring(R, parent(p), cached)) where {T <: RingElement}
   base_ring(parent) != R && error("Base rings do not match.")
   return _map(R, p, parent)
end

@doc raw"""
    change_coefficient_ring(R::Ring, p::MPolyRingElem{<: RingElement}; parent::MPolyRing, cached::Bool=true)

Return the polynomial obtained by coercing the non-zero coefficients of `p`
into `R`.

If the optional `parent` keyword is provided, the polynomial will be an
element of `parent`. The caching of the parent object can be controlled
via the `cached` keyword argument.
"""
function change_coefficient_ring(R::Ring, p::MPolyRingElem{T}; cached::Bool=true, parent::AbstractAlgebra.MPolyRing = _change_mpoly_ring(R, parent(p), cached)) where {T <: RingElement}
  return change_base_ring(R, p, cached = cached, parent = parent)
end

################################################################################
#
#  Map
#
################################################################################

@doc raw"""
    map_coefficients(f, p::MPolyRingElem{<: RingElement}; parent::MPolyRing)

Transform the polynomial `p` by applying `f` on each non-zero coefficient.

If the optional `parent` keyword is provided, the polynomial will be an
element of `parent`. The caching of the parent object can be controlled
via the `cached` keyword argument.
"""
function map_coefficients(f, p::MPolyRingElem; cached = true, parent::AbstractAlgebra.MPolyRing = _change_mpoly_ring(AbstractAlgebra.parent(f(zero(base_ring(p)))), AbstractAlgebra.parent(p), cached))
   return _map(f, p, parent)
end

function _map(g, p::MPolyRingElem, Rx)
   cvzip = zip(coefficients(p), exponent_vectors(p))
   M = Generic.MPolyBuildCtx(Rx)
   for (c, v) in cvzip
      push_term!(M, g(c), v)
   end

   return finish(M)
end

###############################################################################
#
#   Parent object overload (with coercion)
#
###############################################################################

function (S::MPolyRing)(b::Vector, m::Vector{Vector{T}}) where {T}
   R = base_ring(S)
   return S(map(R, b)::Vector{elem_type(R)}, m)
end

###############################################################################
#
#   Random elements
#
###############################################################################

function rand_ordering(rng::AbstractRNG=Random.GLOBAL_RNG)
   i = rand(rng, 1:3)
   if i == 1
      return :lex
   elseif i == 2
      return :deglex
   else
      return :degrevlex
   end
end

RandomExtensions.maketype(S::AbstractAlgebra.MPolyRing, _, _, _) = elem_type(S)

function RandomExtensions.make(S::AbstractAlgebra.MPolyRing, term_range::UnitRange{Int},
                               exp_bound::UnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, exp_bound, vs[1])
   else
      make(S, term_range, exp_bound, make(R, vs...))
   end
end

function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make4{
                 <:RingElement,<:AbstractAlgebra.MPolyRing,UnitRange{Int},UnitRange{Int}}})
   S, term_range, exp_bound, v = sp[][1:end]
   f = S()
   g = gens(S)
   R = base_ring(S)
   for i = 1:rand(rng, term_range)
      term = S(1)
      for j = 1:length(g)
         term *= g[j]^rand(rng, exp_bound)
      end
      term *= rand(rng, v)
      f += term
   end
   return f
end

function rand(rng::AbstractRNG, S::AbstractAlgebra.MPolyRing,
              term_range::UnitRange{Int}, exp_bound::UnitRange{Int}, v...)
   rand(rng, make(S, term_range, exp_bound, v...))
end

function rand(S::AbstractAlgebra.MPolyRing, term_range, exp_bound, v...)
   rand(GLOBAL_RNG, S, term_range, exp_bound, v...)
end

###############################################################################
#
#   polynomial_ring constructor
#
###############################################################################

@doc raw"""
    polynomial_ring(R::Ring, s::Vector{T}; cached::Bool = true, ordering::Symbol = :lex) where T <: VarName

Given a base ring `R` and a vector `s` of variable names $x1, x2, \dots$ specifying
how the generators (variables) should be printed, return a tuple `S, [x1, x2, ...]`
representing the new polynomial ring $S = R[x1, x2, ...]$ and the generators
$x1, x2, \dots$ of the polynomial ring.

Mathematically the object `S` depends only on `R` and `x1, x2, ...` and by
default it will be cached, i.e., if `polynomial_ring` is invoked again with the
same arguments, the same (*identical*) ring is returned. Setting the optional
argument `cached` to `false` ensures a distinct new ring is returned, and will
also prevent it from being cached.

The `ordering` of the polynomial ring can be one of `:lex`, `:deglex` or `:degrevlex`.
"""
polynomial_ring(R::Ring, s::Union{Tuple{Vararg{T}}, AbstractVector{T}}; kw...) where
      T<:VarName =
   polynomial_ring(R, Symbol[Symbol(x) for x in s]; kw...)

function polynomial_ring(R::Ring, s::Vector{Symbol}; kw...)
   S = polynomial_ring_only(R, s; kw...)
   (S, gens(S))
end

@doc raw"""
    polynomial_ring_only(R::Ring, s::Vector{Symbol}; ordering::Symbol=:lex, cached::Bool=true)

Like [`polynomial_ring(R::Ring, s::Vector{Symbol})`](@ref) but return only the
multivariate polynomial ring.
"""
polynomial_ring_only(R::T, s::Vector{Symbol}; ordering::Symbol=:lex, cached::Bool=true) where T<:Ring =
   mpoly_ring_type(T)(R, s, ordering, cached)

# Alternative constructors

@doc raw"""
    polynomial_ring(R::Ring, n::Int, s::VarName = :x; cached, ordering)

Given a symbol, string or character `s` and a number of variables `n` will
do the same as the first constructor except that the variables will be
automatically numbered. For example if `s` is the string `x` and `n = 3` then
the variables will print as `x1`, `x2`, `x3`.
"""
polynomial_ring(R::Ring, n::Int, s::VarName=:x; kw...) =
   polynomial_ring(R, Symbol.(s, 1:n); kw...)

MPolyRing(R::Ring, n::Int) =
   polynomial_ring_only(R, Symbol.(:x, 1:n); cached=false)

################################################################################
#
#  Fancy macro
#
################################################################################

function build_names(prefix, indices...)
   map(i -> "$(prefix)[$(join(i, ","))]", Iterators.product(indices...))
end

function build_variable(arg::Symbol)
   t = gensym()
   return t, :($(esc(t)) = String[$"$arg"])
end

function build_variable(arg::Expr)
   isa(arg, Expr) || error("Expected $arg to be a variable name")
   Base.Meta.isexpr(arg, :ref) || error("Expected $arg to be of the form varname[idxset]")
   (2 ≤ length(arg.args)) || error("Expected $arg to have at least one index set")
   varname = arg.args[1]
   prefix = string(varname)
   t = gensym()
   return t, :($(esc(t)) = build_names($prefix, $(esc.(arg.args[2:end])...)))
end

function build_variables_strings(args)
   names = Symbol[]
   exprs = Expr[]
   for arg in args
     name_var, define_names = build_variable(arg)
     push!(exprs, define_names)
     #push!(exprs, :(print($(esc(name_var)))))
     push!(names, name_var)
   end
   return names, exprs
end

macro polynomial_ring(R, args...)
   names, exprs = build_variables_strings(args)
   all_names = gensym()
   push!(exprs, :($(esc(all_names)) = String[]))
   for t in names
      push!(exprs, :(append!($(esc(all_names)), reshape($(esc(t)), length($(esc(t)))))))
   end
   ring1 = gensym()
   ring2 = gensym()
   push!(exprs, :($(Expr(:tuple, esc(ring1), esc(ring2))) =
                 AbstractAlgebra.polynomial_ring($(esc(R)), $(esc(all_names)))))
   vars = Symbol[]
   k = gensym()
   push!(exprs, :($(esc(k)) = 0))
   for (i, t) in enumerate(names)
      var_sym = gensym()
      if args[i] isa Symbol
         push!(exprs, :($(esc(var_sym)) = ($(esc(ring2)))[$(esc(k)) + 1]))
      else
         push!(exprs, :($(esc(var_sym)) = elem_type($(esc(ring1)))[($(esc(ring2)))[$(esc(k)) + i] for (i,_) in enumerate($(esc(t)))]))
      end
      push!(vars, var_sym)
      push!(exprs, :($(esc(k)) = $(esc(k)) + length($(esc(t)))))
   end
   res = :($(foldl((x,y) -> :($x; $y), exprs, init=:())); $(Expr(:tuple, esc(ring1), esc.(vars)...)))
   return res
end
