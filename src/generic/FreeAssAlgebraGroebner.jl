###############################################################################
#
#   FreeAssAlgebraGroebner.jl : free associative algebra R<x1,...,xn> Groebner basis
#
###############################################################################
#include("../AbstractTypes.jl")
#include("FreeAssAlgebra.jl")
include("FreeAssAhoCorasick.jl")
export groebner_basis, interreduce!

using DataStructures

const groebner_debug_level = 0
const Monomial = Vector{Int}

abstract type Obstruction{T} end 
"""
represents the overlap of the leading term of two polynomials
"""

struct ObstructionTriple{T} <: Obstruction{T} 
        first_poly::FreeAssAlgElem{T}
        second_poly::FreeAssAlgElem{T}
        pre_and_suffixes::NTuple{4, Monomial}
end

function FreeAssAlgElem{T}(R::FreeAssAlgebra{T}, mon::Monomial) where T
        return FreeAssAlgElem{T}(R, [one(base_ring(R))], [mon], 1)
end

"""
takes an obstruction triple (p, q, o) and returns the common multiple 
of the leading terms of p and q defined by o
TODO documentation
"""
function common_multiple_leading_term(ot::ObstructionTriple{T}) where T
    return FreeAssAlgElem{T}(parent(ot.first_poly), ot.pre_and_suffixes[1]) * FreeAssAlgElem{T}(parent(ot.first_poly), _leading_word(ot.first_poly)) * FreeAssAlgElem{T}(parent(ot.first_poly), ot.pre_and_suffixes[2])
end

function s_polynomial(ot::ObstructionTriple{T}) where T
    first_term = FreeAssAlgElem{T}(parent(ot.first_poly), ot.pre_and_suffixes[1]) * ot.first_poly * FreeAssAlgElem{T}(parent(ot.first_poly), ot.pre_and_suffixes[2])
    second_term = FreeAssAlgElem{T}(parent(ot.first_poly), ot.pre_and_suffixes[3]) * ot.second_poly * FreeAssAlgElem{T}(parent(ot.first_poly), ot.pre_and_suffixes[4])
       return inv(leading_coefficient(ot.first_poly)) * first_term - inv(leading_coefficient(ot.second_poly)) * second_term
end

# skip all of the extra length-checking
function _leading_word(a::FreeAssAlgElem{T}) where T
   return a.exps[1]
end

function gb_divides_leftmost_aho_corasick(a::Word, aut::AhoCorasickAutomaton)
    match = search(aut, a)
    if isnothing(match)
        return (false, [], [], -1)
    end
    return (true, a[1:match[1] - length(match[2])], a[match[1] + 1:length(a)], match[2][1])
end

# implementation of the normal form function using aho corasick to check for all groebner basis elements in parallel
function normal_form(
        f::FreeAssAlgElem{T},
        g::Vector{FreeAssAlgElem{T}},
        aut::AhoCorasickAutomaton
    ) where T
    R = parent(f)
    rexps = Vector{Int}[]
    rcoeffs = T[]
    while length(f) > 0
        ok, left, right, match_index = gb_divides_leftmost_aho_corasick(f.exps[1], aut)
        if ok
            qi = divexact(f.coeffs[1], g[match_index].coeffs[1])
            f = _sub_rest(f, mul_term(qi, left, g[match_index], right), 1)
        else
            push!(rcoeffs, f.coeffs[1])
            push!(rexps, f.exps[1])
            f = FreeAssAlgElem{T}(R, f.coeffs[2:end], f.exps[2:end], length(f)-1)
        end
    end
    return FreeAssAlgElem{T}(R, rcoeffs, rexps, length(rcoeffs))
end


# normal form with leftmost word divisions
function normal_form(
   f::FreeAssAlgElem{T},
   g::Vector{FreeAssAlgElem{T}},
   suffix_match_vectors::Vector{Vector{Int}}
) where T <: FieldElement
   R = parent(f)
   s = length(g)
   rcoeffs = T[]
   rexps = Vector{Int}[]
   while length(f) > 0
      i = 1
   @label again
        ok, ml, mr = word_divides_leftmost(f.exps[1], g[i].exps[1], suffix_match_vectors[i])
        if !ok && i < s
         i += 1
         @goto again
      end
      if ok
         qi = divexact(f.coeffs[1], g[i].coeffs[1])
         f = _sub_rest(f, mul_term(qi, ml, g[i], mr), 1) # enforce lt cancelation
      else
         push!(rcoeffs, f.coeffs[1])
         push!(rexps, f.exps[1])
         f = FreeAssAlgElem{T}(R, f.coeffs[2:end], f.exps[2:end], length(f)-1)
      end
   end
   r = FreeAssAlgElem{T}(R, rcoeffs, rexps, length(rcoeffs))
   return r
end

# weak normal form with leftmost word divisions
function normal_form_weak(
   f::FreeAssAlgElem{T},
   g::Vector{FreeAssAlgElem{T}}
) where T <: FieldElement
   R = parent(f)
   s = length(g)
   while length(f) > 0
      i = 1
   @label again
      ok, ml, mr = word_divides_leftmost(f.exps[1], g[i].exps[1])
      if !ok && i < s
         i += 1
         @goto again
      end
      if ok
         qi = divexact(f.coeffs[1], g[i].coeffs[1])
         f = _sub_rest(f, mul_term(qi, ml, g[i], mr), 1) # enforce lt cancelation
      else
         break
      end
   end
   return f
end

function interreduce!(g::Vector{FreeAssAlgElem{T}}) where T
   i = 1
   while length(g) > 1 && length(g) >= i
      r = normal_form(g[i], g[1:end .!= i])
      if iszero(r)
         deleteat!(g, i)
      elseif g[i] != r
         g[i] = r
         i = 1
      else
         i += 1
      end
   end
   return g
end

## checks whether there is an overlap between a and b at position i of b
#  such that b[i:length(b)] = a[1:length(b)-i]
function check_left_overlap(a::Vector{Int}, b::Vector{Int}, i::Int)
   for j in 0:length(b)-i
      if j >= length(a)
         return false # this is a center overlap
      end
      if b[i+j] != a[j+1]
         return false
      end
   end
   return true
end

###
# find all non-trivial left-obstructions of a and b
# i.e. all words w_1 and w_2^' s.t. w_1 a = b w_2^'
# where length(w_1) < length(b) and length(w_2^') < length(a)
# the return vector is of the form [(w_1, w_2^'), ...]
# if w_1 or w_2^' is empty, the corresponding obstruction is not returned
function left_obstructions(a::Vector{Int}, b::Vector{Int})
   v = Tuple{Vector{Int}, Vector{Int}}[];
   for i in 2:length(b)
      if check_left_overlap(a, b, i)
         if length(b)-i + 2 <= length(a) # w_2^' should not be empty!
            push!(v, (b[1:i-1], a[length(b)-i + 2:length(a)]))
         end
      end
   end
   return v
end


###
# find all non-trivial right-obstructions of a and b
# i.e. all words w_2 and w_1^' s.t. a w_1^' = w_2 b
# where length(w_1^') < length(b) and length(w_2) < length(a)
# the return vector is of the form [(w_2, w_1^'), ...]
# if w_1^' or w_2 is empty, the corresponding obstruction is not returned
function right_obstructions(a::Vector{Int}, b::Vector{Int})
   return left_obstructions(b, a)
end

###
# check, whether a is a true subword of b at index i
function check_center_overlap(a::Vector{Int}, b::Vector{Int}, i::Int)
   for j in 1:length(a)
      if i+j -1 > length(b)
         return false
      end
      if a[j] != b[i + j - 1]
         return false
      end
   end
   return true
end

function center_obstructions_first_in_second(a::Vector{Int}, b::Vector{Int})
   v = Tuple{Vector{Int}, Vector{Int}}[]
   for i in 1:length(b)
      if check_center_overlap(a, b, i)
         push!(v, (b[1:i-1], b[i + length(a): length(b)]))
      end
   end
   return v
end

##
# return all center obstructions of a and b, i.e. all (w_i, w_i^')
# such that either
# w_i a w_i^' = b
# or
# w_i b w_i^' = a
# either or both of w_i and w_i^' can be empty
function center_obstructions(a::Vector{Int}, b::Vector{Int})
   if length(a) > length(b)
      return center_obstructions_first_in_second(b, a)
   else
      return center_obstructions_first_in_second(a, b)
   end
end

# all non-trivial ones
function obstructions(a::Vector{Int}, b::Vector{Int})
   one = Int[] # the empty word
   res = Tuple{Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}}[]
   for x in center_obstructions_first_in_second(b, a)
      push!(res, (one, one, x[1], x[2]))
   end
   for x in center_obstructions_first_in_second(a, b)
      push!(res, (x[1], x[2], one, one))
   end
   for x in left_obstructions(a, b)
      push!(res, (x[1], one, one, x[2]))
   end
   for x in left_obstructions(b, a)
      push!(res, (one, x[2], x[1], one))
   end
   for x in res
      @assert vcat(x[1], a, x[2]) == vcat(x[3], b, x[4])
   end
   return res
end

# all non-trivial self obstructions
function obstructions(a::Vector{Int})
   one = Int[] # the empty word
   res = Tuple{Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}}[]
   for x in left_obstructions(a, a)
      push!(res, (one, x[2], x[1], one))
   end
   for x in res
      @assert vcat(x[1], a, x[2]) == vcat(x[3], a, x[4])
   end
   return res
end


# check whether w_2 = v w_1 for some word v
function is_subword_right(w_1::Vector{Int}, w_2::Vector{Int})
   if length(w_1) > length(w_2)
      return false
   end
   for i in 0:length(w_1) - 1
      if w_1[length(w_1) - i] != w_2[length(w_2) - i]
         return false
      end
   end
   return true
end


# check whether w_2 = w_1 v for some word v
function is_subword_left(w_1::Vector{Int}, w_2::Vector{Int})
   if length(w_1) > length(w_2)
      return false
   end
   for i in 1:length(w_1)
      if w_1[i] != w_2[i]
         return false
      end
   end
   return true
end


###
# check if for obs1 = (w_i, w_i^'; u_j, u_j^') and obs2 = (w_k, w_k^'; v_l, v_l^')
# it holds that u_j == w v_l and u_j^' = v_l^' w^' for some w, w^'
# i.e. if obs2 is a subobstruction of obs1
# both w and w^' might be empty
function is_subobstruction(obs1::NTuple{4, Vector{Int}}, obs2::NTuple{4, Vector{Int}})
#        if length(obs2[3]) > length(obs1[3]) || length(obs2[4]) > length(obs1[4])
#                return false
#        end
   if is_subword_right(obs2[3], obs1[3]) && is_subword_left(obs2[4], obs1[4])
      return true
   else
      return false
   end
end

# check whether there exists a w^'' such that
# w1 LM(g1) w2 = w1 LM(g1) w^'' LM(g2) u2
function has_overlap(g1, g2, w1, w2, u1, u2)
   lw1 = _leading_word(g1)
   lw2 = _leading_word(g2)
   concatenated_word = vcat(w1, lw1, w2)
   for i in 1:length(w1)
      c = popfirst!(concatenated_word)
      @assert c == w1[i]
   end
   for i in 1:length(lw1)
      c = popfirst!(concatenated_word)
      @assert c == lw1[i]
   end
   for j in 0:length(u2)-1
      c = pop!(concatenated_word)
      @assert c = u2[length(u2)-j]
   end
   if length(concatenated_word) < length(lw2)
      return false
   end
   return is_subword_right(lw2, concatenated_word) # TODO maybe just comparing lengths should be sufficient
end

function is_redundant(# TODO do we need g in the signature?
   obs::NTuple{4, Vector{Int}},
   obs_index::Int,
   s::Int,
   B::Matrix{Vector{NTuple{4, Vector{Int}}}},
   g::Vector{FreeAssAlgElem{T}}
) where T
   # cases 4b + 4c
   for j in 1:size(B)[1]
      for k in 1:length(B[j, s])
         if is_subobstruction(obs, B[j, s][k])
            # case 4b
            if length(obs[3]) - length(B[j, s][k][3]) + length(obs[4]) - length(B[j,s][k][4]) > 0
               return true
            # case 4c
            elseif obs_index > j
               return true
            elseif obs_index == j &&
                   length(obs[3]) - length(B[j, s][k][3]) + length(obs[4]) - length(B[j,s][k][4]) == 0 &&
                   word_gt(obs[1], B[j, s][k][1])
               return true
            end
         end
      end
   end
   # case 4d
   # size(B) should be (s, s)
   # we want to iterate over all B[i, obs_index] with i < s
   for i in 1:size(B)[1]-1
      for k in 1:length(B[i, obs_index])
         if is_subword_right(B[i, obs_index][k][3], obs[1]) && is_subword_left(B[i, obs_index][k][4], obs[2])
            if groebner_debug_level > 0
               show(obs)
               show(B[i, obs_index][k])
            end
            u1 = copy(obs[1])
            u2 = copy(obs[2])
            v1 = copy(B[i, obs_index][k][3])
            v2 = copy(B[i, obs_index][k][4])
            for i in 1:length(v1)
               pop!(u1)
            end
            for j in 1:length(v2)
               popfirst!(u2)
            end
            @assert word_cmp(vcat(u1, v1), obs[1]) == 0
            @assert word_cmp(vcat(v2, u2), obs[2]) == 0
            if !has_overlap(g[i], g[s], vcat(u1, B[i, obs_index][k][1]), vcat(B[i, obs_index][k][2], u2), obs[3], obs[4])
               return true
            end
         end
      end
   end
   return false
end

## s is the index of the newly added layer of obstructions in B
function remove_redundancies!(
   B::Matrix{Vector{NTuple{4, Vector{Int}}}},
   s::Int,
   g::Vector{FreeAssAlgElem{T}}
) where T
   del_counter = 0
   for i in 1:s
      k = 1
      while k <= length(B[i, s])
         if is_redundant(B[i, s][k], i, s, B, g)
            deleteat!(B[i, s], k)
            del_counter += 1
         else
            k += 1
         end
      end
   end
   # TODO case 4e from Thm 4.1 in Kreuzer Xiu
end

function get_obstructions(g::Vector{FreeAssAlgElem{T}}) where T
    s = length(g)
    result = PriorityQueue{Obstruction{T}, FreeAssAlgElem{T}}()
    for i in 1:s, j in 1:i
        if i == j
            obs = obstructions(_leading_word(g[i]))
        else
            obs = obstructions(_leading_word(g[i]), _leading_word(g[j]))
        end
        for o in obs
            triple = ObstructionTriple{T}(g[i], g[j], o)
            enqueue!(result, triple, common_multiple_leading_term(triple))
        end
    end
    # TODO maybe here some redundancies can be removed too, check Kreuzer Xiu
    if groebner_debug_level > 0
        obstr_count = length(result)
        println("$obstr_count many obstructions")
    end
    return result
end


function add_obstructions!(
    obstruction_queue::PriorityQueue{Obstruction{T}, FreeAssAlgElem{T}},
    g::Vector{FreeAssAlgElem{T}}
) where T
    s = length(g)
    for i in 1:s
        if i == s
            obs = obstructions(_leading_word(g[i]))
        else
            obs = obstructions(_leading_word(g[i]), _leading_word(g[s]))
        end
        for o in obs
            triple = ObstructionTriple{T}(g[i], g[s], o)
            enqueue!(obstruction_queue, triple, common_multiple_leading_term(triple))
        end
    end
    #remove_redundancies!(new_B, s, g) TODO match remove_redundancies to new types
end


function groebner_basis_buchberger(
   g::Vector{FreeAssAlgElem{T}},
   reduction_bound = typemax(Int)::Int
) where T <: FieldElement

   g = copy(g)
   checked_obstructions = 0
   nonzero_reductions = 0
   # compute the aho corasick automaton
   # to make normal form computation more efficient
   aut = AhoCorasickAutomaton([g_i.exps[1] for g_i in g])

   # step 1
   obstruction_queue = get_obstructions(g) 
   while !isempty(obstruction_queue)
      obstruction = dequeue!(obstruction_queue)
      # step3 
      S = s_polynomial(obstruction)
      Sp = normal_form(S, g, aut) # or normal_form_weak
      if groebner_debug_level > 0
          checked_obstructions += 1
          if checked_obstructions % 5000 == 0
            println("checked $checked_obstructions obstructions")
         end
      end
      if iszero(Sp)
         continue
      end
      nonzero_reductions += 1
      # step4
      push!(g, Sp)
      aut = AhoCorasickAutomaton([g_i.exps[1] for g_i in g])

      if groebner_debug_level > 0
         println("adding new obstructions! checked $checked_obstructions so far")
      end
      if nonzero_reductions >= reduction_bound
              return g
      end
      add_obstructions!(obstruction_queue, g)
   end
   return g
end

function groebner_basis(
   g::Vector{FreeAssAlgElem{T}},
   reduction_bound = typemax(Int)::Int
) where T <: FieldElement
   return groebner_basis_buchberger(g, reduction_bound)
end

