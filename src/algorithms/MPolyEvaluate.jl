###############################################################################
#
#   MPolyEvaluate.jl : Various algorithms for evaluating multivariates
#
###############################################################################

function add!(z::T, a, b)::T where T
  return a + b
end

function mul!(z::T, a, b)::T where T
  return a*b
end

function pow!(z::T, a, b)::T where T
  return a^b
end

# a *= x^p, using ta and tx as temps. Must return the modified a, ta, tx
function mulpow!(a::T, x::S, p, ta::T, tx::S) where {T, S}
  tx = pow!(tx, x, p)
  ta = mul!(ta, a, tx)
  (a, ta) = (ta, a)
  return (a, ta, tx)
end

function mulpow!(a::T, x::S, p, ta::T, tx::S) where {T <: MPolyRingElem, S <: MPolyRingElem}
  if p > 1 && (length(x) < 2 || length(a) < p*length(x))
    tx = pow!(tx, x, p)
    ta = mul!(ta, a, tx)
    (a, ta) = (ta, a)
  else
    while p > 0
      ta = mul!(ta, a, x)
      (a, ta) = (ta, a)
      p -= 1
    end
  end
  return (a, ta, tx)
end

###############################################################################
#
#   evaluate_horner_lex
#
###############################################################################

# The usual univariate horner rule. Since Bexps is in lex order,
# recursion is easy, and the depth is bounded by nvars
function _horner_lex_rec(
  res::Vector,
  ctxA::Ring,
  Bcoeffs, Bexps::Vector{Vector{Int}},
  Bstart::Int, Bstop::Int, # ignore all terms outside this interval [Bstart, Bstop]
  var::Int, # ignore all variables of index < var
  C::Vector,
  ctxC::Ring,
  ta, tc) # temps

  @assert Bstart <= Bstop
  len = Bstop + 1 - Bstart

  if len < 2
    # evaluate a monomial
    first = true
    for v in var:length(Bexps[Bstart])
      Bexps[Bstart][v] == 0 && continue
      if first
        tc = pow!(tc, C[v], Bexps[Bstart][v])
        res[var] = mul!(res[var], Bcoeffs[Bstart], tc)
      else
        res[var], ta, tc = mulpow!(res[var], C[v], Bexps[Bstart][v], ta, tc)
      end
      first = false
    end
    if first
      res[var] = deepcopy(ctxA(Bcoeffs[Bstart]))  # res[.] should be mutable
      first = false
    end
    return (ta, tc)
  end

  i = Bstart
  e = Bexps[i][var]
  res[var] = zero!(res[var])

@label next_term

  i += 1
  if i > Bstop
    nexte = 0
  elseif Bexps[i][var] == e
    @goto next_term
  else
    nexte = Bexps[i][var]
  end

  (ta, tc) = _horner_lex_rec(res, ctxA, Bcoeffs, Bexps, Bstart, i-1, var+1, C, ctxC, ta, tc)
  res[var] = addeq!(res[var], res[var+1])
  res[var], ta, tc = mulpow!(res[var], C[var], e - nexte, ta, tc)

  if i > Bstop
    return (ta, tc)
  end

  Bstart = i
  e = nexte
  @goto next_term
end

function evaluate_horner_lex(B::MPolyRingElem, C::Vector{<:RingElement})
  @assert nvars(parent(B)) <= length(C)
  ctxC = parent(C[1])
  Bcoeffs = collect(coefficients(B))
  ctxA = parent(one(coefficient_ring(B)) * one(ctxC))
  Bexps = collect(exponent_vectors(B))
  n = length(Bexps)
  if n < 1
    return zero(ctxA)
  end
  p = sortperm(Bexps, rev=true)
  if p != 1:n
    Bcoeffs = [Bcoeffs[p[i]] for i in 1:n]
    Bexps   = [Bexps[p[i]] for i in 1:n]
  end
  res = elem_type(ctxA)[ctxA() for i in 1:nvars(parent(B))+2]
  _horner_lex_rec(res, ctxA, Bcoeffs, Bexps, 1, n, 1, C, ctxC, ctxA(), ctxC())
  return res[1]
end

###############################################################################
#
#   evaluate_horner
#
###############################################################################

#=
This conversion to Horner form can be stated as recursive. However, the call
stack has depth proportional to the length of the input polynomial in the worst
case. Therefore, the following iterative spaghetti ensues.

HornerForm(f):
  if f is simple to evaluate (choice in else branch doesn't exist)
    return eval(f)
  else
    choose a variable v and the smallest nonzero exponent e appearing in f
    write f = q * v^e + r  where r is independent of the variable v
    return  HornerForm(q) * v^e + HornerForm(r)

TODO: change all zero-based indices into natural one-based indices

Since the evaluation proceeds in a top-down fashion (instead of bottom-up)
the overhead of generating all of these q's and r's is quadratic in length(f)
in the worst case.
=#
struct stack_entry
  f::Int
  r::Int
  v_var::Int
  v_exp::Int
  ret::Int
end

# obtain a count of the number of terms containing each variable
function _term_counts!(
  counts::Vector{Int},  # modified
  mdegs::Vector{Int},   # modified
  nvars::Int,
  f::Int,
  Bexps::Vector{Vector{Int}},
  Blist::Vector{Int}
)
  for i in 0:nvars-1
    counts[1+i] = 0
    mdegs[1+i] = -1
  end
  j = f
  while j != -1
    for i in 0:nvars-1
      if Bexps[1+j][1+i] != 0
        counts[1+i] += 1
        if mdegs[1+i] < 0 || mdegs[1+i] > Bexps[1+j][1+i]
          mdegs[1+i] = Bexps[1+j][1+i]
        end
      end
    end
    j = Blist[1+j]
  end

  totalcounts = 0
  maxcounts = 0
  v = -1
  for i in 0:nvars-1
    maxcounts = max(maxcounts, counts[1+i])
    totalcounts += counts[1+i]
    if counts[1+i] != 0
      v = i
    end
  end

  return totalcounts, maxcounts, v
end

# split into q and r
function _split_qr!(
  counts::Vector{Int},
  mdegs::Vector{Int},
  nvars::Int,
  f::Int,
  Bexps::Vector{Vector{Int}}, # modified
  Blist::Vector{Int},         # modified
  maxcounts::Int
)
  # pick best power to pull out
  k = 0
  if maxcounts == 1
    score = -1
    for i in 0:nvars-1
      if counts[1+i] == 1 && (score < 0 || mdegs[1+i] < score)
        @assert mdegs[1+i] > 0
        score = mdegs[1+i]
        k = i
      end
    end
  else
    score = 0
    for i in 0:nvars-1
      if counts[1+i] > 1
        @assert mdegs[1+i] > 0
        if mdegs[1+i]*(counts[1+i]-1) > score
          score = mdegs[1+i]*(counts[1+i]-1)
          k = i
        end
      end
    end
  end

  # scan f and split into q and v with f = q*v + r
  r = -1
  cur = f
  f_prev = -1
  r_prev = -1
  while cur != -1
    next = Blist[1+cur]
    if Bexps[1+cur][1+k] == 0
      if f_prev == -1
        f = Blist[1+cur]
      else
        Blist[1+f_prev] = Blist[1+cur]
      end
      if r_prev == -1
        r = cur
      else
        Blist[1+r_prev] = cur
      end
      Blist[1+cur] = -1
      r_prev = cur
    else
      # mdegs[k] should be minimum non zero exponent
      @assert Bexps[1+cur][1+k] >= mdegs[1+k]
      Bexps[1+cur][1+k] -= mdegs[1+k]
      f_prev = cur
    end
    cur = next
  end
  return f, r, k
end

# evaluate B at vars = C
# arithmetic between ctxB and ctxC should be in ctxA
function _evaluate_horner_non_rec(
  ctxA::Ring,
  Bcoeffs,
  Bexps::Vector{Vector{Int}},
  C::Vector,
  ctxC::Ring
)

  Blen = length(Bcoeffs)
  if Blen < 1
    return zero(ctxA)
  end

  nvars = length(Bexps[1])
  if nvars < 1
    return ctxA(Bcoeffs[1+0])
  end

  counts = zeros(Int, nvars)
  mdegs = zeros(Int, nvars)
  stack = Vector{stack_entry}()

  # regs and rtypes: Hold the return values of HornerForm.
  #   rtypes is used like a stack with stack pointer rp.
  #   If rtypes[1+rp] == -1, then the return is in regs[1+rp]
  #   Otherwise, the return is in Bcoeffs[1+rtypes[1+rp]]
  rp = 0
  rtypes = zeros(Int, nvars + 1)
  regs = elem_type(ctxA)[ctxA() for i in 1:nvars]
  tempC = ctxC()
  tempA = ctxA()

  # the polynomials will be stored as null-terminated linked lists, null = -1
  Blist = zeros(Int, Blen)
  for i in 0:Blen-2
    Blist[1+i] = i+1
  end
  Blist[1+Blen-1] = -1

  # start with f = B
  push!(stack, stack_entry(0, 0, 0, 0, 0))

@label horner_form  # beginning of HornerForm(f)

  f = stack[end].f
  @assert f != -1

  totalcounts, maxcounts, v = _term_counts!(counts, mdegs, nvars, f, Bexps, Blist)

  # handle simple cases
  if totalcounts == 0
    @assert Blist[1+f] == -1  # f should have only one term
    rtypes[1+rp] = f
    @goto horner_form_return
  elseif totalcounts == 1
    @assert Bexps[1+f][1+v] > 0
    if elem_type(ctxA) == elem_type(ctxC)
      regs[1+rp] = pow!(regs[1+rp], C[1+v], Bexps[1+f][1+v])
      regs[1+rp] = mul!(regs[1+rp], regs[1+rp], Bcoeffs[1+f])
    else
      tempC = pow!(tempC, C[1+v], Bexps[1+f][1+v])
      regs[1+rp] = mul!(regs[1+rp], tempC, Bcoeffs[1+f])
    end
    if Blist[1+f] != -1       # f has a second term that should be a scalar
      @assert Bexps[1+Blist[1+f]][1+v] == 0
      regs[1+rp] = add!(regs[1+rp], regs[1+rp], Bcoeffs[1+Blist[1+f]])
    end
    rtypes[1+rp] = -1
    @goto horner_form_return
  end

  f, r, k = _split_qr!(counts, mdegs, nvars, f, Bexps, Blist, maxcounts)

  # call frame for the remainder
  stack[end] = stack_entry(stack[end].f, r, k, mdegs[1+k], stack[end].ret)

  # convert the quotient
  push!(stack, stack_entry(f, 0, 0, 0, 1))
  @goto horner_form
@label horner_form_ret_1

  # convert the remainder
  r = stack[end].r

  if r != -1
    rp += 1
    @assert 0 <= rp <= nvars
    push!(stack, stack_entry(r, 0, 0, 0, 2))
    @goto horner_form
@label horner_form_ret_2

    if rtypes[1+rp-1] == -1 && rtypes[1+rp] == -1
      # both quotient and remainder are polynomials
      regs[1+rp-1], tempA, tempC = mulpow!(regs[1+rp-1], C[1+stack[end].v_var], stack[end].v_exp, tempA, tempC)
      tempA = add!(tempA, regs[1+rp-1], regs[1+rp])
      regs[1+rp-1], tempA = tempA, regs[1+rp-1]
    elseif rtypes[1+rp-1] == -1 && rtypes[1+rp] != -1
      # quotient is a polynomial, remainder is a scalar
      regs[1+rp-1], tempA, tempC = mulpow!(regs[1+rp-1], C[1+stack[end].v_var], stack[end].v_exp, tempA, tempC)
      regs[1+rp-1] = add!(regs[1+rp-1], regs[1+rp-1], Bcoeffs[1+rtypes[1+rp]])
    elseif rtypes[1+rp-1] != -1 && rtypes[1+rp] == -1
      # quotient is a scalar, remainder is a polynomial
      tempC = pow!(tempC, C[1+stack[end].v_var], stack[end].v_exp)
      tempA = mul!(tempA, tempC, Bcoeffs[1+rtypes[1+rp-1]])
      regs[1+rp-1] = add!(regs[1+rp-1], regs[1+rp], tempA)
    else
      # done by simple case
      @assert false
    end
    rp -= 1
    @assert 0 <= rp <= nvars
  else
    # remainder is zero
    @assert rtypes[1+rp] == -1
    # quotient is a polynomial
    regs[1+rp], tempA, tempC = mulpow!(regs[1+rp], C[1+stack[end].v_var], stack[end].v_exp, tempA, tempC)
  end
  rtypes[1+rp] = -1

@label horner_form_return
  # two possible recursive call sites
  ret = stack[end].ret
  pop!(stack)
  if ret == 1
    @goto horner_form_ret_1
  elseif ret == 2
    @goto horner_form_ret_2
  end

  @assert rp == 0
  @assert length(stack) == 0

  if rtypes[1+rp] == -1
    return regs[1+rp]
  else
    return ctxA(Bcoeffs[1+rtypes[1+rp]])
  end
end

function evaluate_horner(B::MPolyRingElem, C::Vector{<:RingElement})
  @assert nvars(parent(B)) <= length(C)
  ctxC = parent(C[1])
  Bcoeffs = collect(coefficients(B))
  ctxA = parent(one(coefficient_ring(B)) * one(ctxC))
  Bexps = collect(exponent_vectors(B))::Vector{Vector{Int}}
  return _evaluate_horner_non_rec(ctxA, Bcoeffs, Bexps, C, ctxC)
end


###############################################################################
#
#   evaluate_log
#
###############################################################################

########################################################################
#
#  Evaluation with multivariate logarithmic caching.
#
#  The code does the following. Given a multivariate polynomial
#  f = ∑ₐ cₐ ⋅ xᵃ in n variables and n elements p₁,…,pₙ on which
#  f will be evaluated, iterate through the terms cₐ ⋅ xᵃ of f.
#
#  On first stage, we have a stack with precomputed powers
#  from previous terms where every entry is divisible by the
#  previous one. Given xᵃ, pop from that stack until the top
#  element q = pᵇ divides pᵃ. Proceed with the computation
#  of pᵃ:q.
#
#  Write the exponent a - b as a sum of integer vectors
#  a - b = c⁽¹⁾ + c⁽²⁾+ ... + c⁽ᵐ⁾ where each of the c⁽ʲ⁾ has integer
#  entries (k₁, k₂, …, kₙ) that are either zero or a pure power
#  of 2.
#
#  Then look up each one of the products p₁ᵏ¹ ⋅ p₂ᵏ² ⋅ … ⋅ pₙᵏⁿ
#  in a cache dictionary. If it's not there, these can be easily
#  computed as products of squares of other precomputed entries.
#
#  Once all these lookups have been carried out, assemble the
#  final product for the monomial q ⋅ pᵃ ⁻ ᵇ using the obvious
#  m+1 multiplications and add them to the return list for
#  summation.
#
#  NOTE: In order for the above algorithm to be efficient,
#  one must iterate through the terms of f starting with
#  those of lowest degree first and going higher. Adjust this
#  according to the iterators!
#
#  Moreover, we assume the iteration to proceed with consecutive
#  terms "close to one another" regarding their exponent vectors.
#
#  This is, for instance, the case when the terms are returned in
#  lexicographical order.
#
########################################################################

# change _lowest_bit to _highest_bit in eval_mon to iterate in different order
#function _highest_bit(k::Int)
#  return 1 << (8*sizeof(Int)-1 - leading_zeros(k))
#end

function _lowest_bit(k::Int)
  return k & -k
end

# for evaluating prod_i p[i]^e[i] for a given e and fixed p
mutable struct PowerCache{T}
  p::Vector{T}

  # holds prod_i p[i]^e[i] where each e[i] is zero or a power of 2
  power_cache::Dict{Vector{Int}, T}

  # monomials earlier in the cache divide later ones
  # the sentinel monomial (0,0,...,0) is always present
  cache::Vector{Tuple{Vector{Int}, T}}

  function PowerCache{T}(p) where T
    n = length(p)
    return new{T}(p,
                  Dict([[i==j ? 1 : 0 for i in 1:n] => p[j] for j in 1:n]),
                  Tuple{Vector{Int}, T}[(zeros(Int, n), one(p[1]))]
                 )
  end
end

# should allow e to be mutated !!
function look_up(PC::PowerCache, e::Vector{Int})
  all(x->(x==0), e) && return one(PC.p[1])
  haskey(PC.power_cache, e) && return PC.power_cache[e]
  e = deepcopy(e)   # !!
  if all(x->(x<=1), e)
    power = prod([PC.p[i] for i in 1:length(PC.p) if e[i] == 1])
    PC.power_cache[e] = power
    return power
  end
  c = [k >> 1 for k in e]
  d = e - c
  power = look_up(PC, c)*look_up(PC, d)
  PC.power_cache[e] = power
  return power
end


function eval_mon(PC::PowerCache, a::Vector{Int})
  n = length(PC.p)
  # pop the cache until the last monomial divides a
  while begin b, q = last(PC.cache); e = a - b; any(x->x<0, e) end
    pop!(PC.cache)
  end

  # a = b + e
  if haskey(PC.power_cache, e)
    qq = q*PC.power_cache[e]
    push!(PC.cache, (a, qq))
    return qq
  end

  # for the monomial e = [7,5], use e = (4,4) + (2,0) + (1,1)
  mask = e[1]
  for i in 2:n
    mask |= e[i]
  end
  newe = zeros(Int, n)
  result = q
  while mask > 0
    bit = _lowest_bit(mask)
    mask = mask & ~bit
    for i in 1:n
      newe[i] = e[i] & bit
    end
    result *= look_up(PC, newe)
    b += newe
    push!(PC.cache, (b, result))
  end
  return result
end

function evaluate_log(f::MPolyRingElem, p::Vector{T}; power_cache = PowerCache{T}(p)) where {T}
  n = length(p)
  R = parent(f)
  n == nvars(R) || error("number of components must equal the number of variables")
  kk = coefficient_ring(R)
  S = parent(one(kk)*one(p[1]))

  r = elem_type(S)[]
  i = UInt(0)
  for (c, e) in reverse(collect(zip(coefficients(f), exponent_vectors(f))))
    push!(r, c*eval_mon(power_cache, e))
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
  return isempty(r) ? zero(S) : r[1]
end

# For a list of polynomials f to be evaluated on the same elements p,
# we can pass on the caching dictionary.
function evaluate_log(f::Vector{PolyType}, p::Vector{T}) where {PolyType<:MPolyRingElem, T}
  D = PowerCache{T}(p)
  return [evaluate_log(g, p; power_cache = D) for g in f]
end

