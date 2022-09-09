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

# a *= x^p, using t as temp. Must return the modifed a and t
function mulpow!(a, x, p, t)
  t = pow!(t, x, p)
  a = mul!(a, a, t)
  return (a, t)
end


#=
The conversion to Horner form can be stated as recursive. However, the call
stack has depth proportional to the length of the input polynomial in the worst
case. Therefore, the following iterative spaghetti ensues.

HornerForm(f):
  if f is simple to evaluate (choice in else branch doesn't exist)
    return eval(f)
  else
    choose a variable v and the smallest nonzero exponent e appearing in f
    write f = q * v^e + r  where r is independent of the variable v
    return  HornerForm(q) * v^e + HornerForm(r)

TODO: change all zero-based indicies into natural one-based indicies

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
function _evaluate_horner_non_rec(
  ctxA::Ring,
  Bcoeffs, Bexps::Vector{Vector{Int}},
  C::Vector,
  ctxC::Ring
)

  Blen = length(Bcoeffs)
  if Blen < 1
    return zero(ctxC)
  end

  nvars = length(Bexps[1])
  if nvars < 1
    return ctxC(Bcoeffs[1+0])
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
  regs = elem_type(ctxC)[ctxC() for i in 1:nvars]
  temp = ctxC()

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
    regs[1+rp] = pow!(regs[1+rp], C[1+v], Bexps[1+f][1+v])
    regs[1+rp] = mul!(regs[1+rp], regs[1+rp], Bcoeffs[1+f])
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
      regs[1+rp-1], temp = mulpow!(regs[1+rp-1], C[1+stack[end].v_var], stack[end].v_exp, temp)
      temp = add!(temp, regs[1+rp-1], regs[1+rp])
      regs[1+rp-1], temp = temp, regs[1+rp-1]
    elseif rtypes[1+rp-1] == -1 && rtypes[1+rp] != -1
      # quotient is a polynomial, remainder is a scalar
      regs[1+rp-1], temp = mulpow!(regs[1+rp-1], C[1+stack[end].v_var], stack[end].v_exp, temp)
      regs[1+rp-1] = add!(regs[1+rp-1], regs[1+rp-1], Bcoeffs[1+rtypes[1+rp]])
    elseif rtypes[1+rp-1] != -1 && rtypes[1+rp] == -1
      # quotient is a scalar, remainder is a polynomial
      temp = pow!(temp, C[1+stack[end].v_var], stack[end].v_exp)
      temp = mul!(temp, temp, Bcoeffs[1+rtypes[1+rp-1]])
      regs[1+rp-1] = add!(regs[1+rp-1], temp, regs[1+rp])
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
    regs[1+rp], temp = mulpow!(regs[1+rp], C[1+stack[end].v_var], stack[end].v_exp, temp)
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
    return ctxC(Bcoeffs[1+rtypes[1+rp]])
  end
end

function evaluate_horner(B::MPolyElem, C::Vector{<:RingElement})
  @assert nvars(parent(B)) <= length(C)
  ctxC = parent(C[1])
  Bcoeffs = collect(coefficients(B))
  ctxA = parent(one(coefficient_ring(B)) * one(ctxC))
  Bexps = collect(exponent_vectors(B))
  return _evaluate_horner_non_rec(ctxA, Bcoeffs, Bexps, C, ctxC)
end

