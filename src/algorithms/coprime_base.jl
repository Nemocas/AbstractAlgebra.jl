################################################################################
#
#  Coprime bases
#
################################################################################

function augment_coprime_base(S::Vector{E}, a::E, start::Int = 1) where E
  i = start
  if is_unit(a)
    return S
  end

  g = a
  new = true

  while i<=length(S) && !isone(a)
    if new
      g = gcd(S[i], a)
      new = false
    else
      g = gcd!(g, S[i], a)
    end
    if is_unit(g)
      i += 1
      continue
    end
    si = divexact(S[i], g)
    a = divexact(a, g)
    if is_unit(si) # g = S[i] and S[i] | a
      continue
    end
    S[i] = si
    if is_unit(a) # g = a and a | S[i]
      a = copy(g)
      continue
    end
    augment_coprime_base(S, copy(g), i)
    continue
  end
  if !is_unit(a)
    push!(S, a)
  end

  return S
end

@doc raw"""
    coprime_base(S::Vector{RingElement}) -> Vector{RingElement}

Returns a coprime base for $S$, i.e. the resulting array contains pairwise coprime objects that multiplicatively generate the same set as the input array.
"""
function coprime_base(S::Vector{E}) where {E <: RingElement}
  return coprime_base_steel(S)
end

function coprime_base_steel(S::Vector)
  @assert !isempty(S)
  T = Array{E}(undef, 1)
  T[1] = S[1]
  for i=2:length(S)
    augment_coprime_base(T, S[i])
  end
  return T
end

@doc raw"""
    coprime_base_push!(S::Vector{RingElem}, a::RingElem) -> Vector{RingElem}

Given an array $S$ of coprime elements, insert a new element, that is, find a
coprime base for `push(S, a)`.
"""
coprime_base_push!(S, a) = augment_coprime_base(S, a)
