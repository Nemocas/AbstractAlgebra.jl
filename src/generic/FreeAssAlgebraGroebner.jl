###############################################################################
#
#   FreeAssAlgebraGroebner.jl : free associative algebra R<x1,...,xn> Groeber basis
#
###############################################################################

export groebner_basis, interreduce!

const groebner_debug_level = 0

# skip all of the extra length-checking
function _leading_word(a::FreeAssAlgElem{T}) where T
   return a.exps[1]
end

# normal form with leftmost word divisions
function normal_form(
   f::FreeAssAlgElem{T},
   g::Vector{FreeAssAlgElem{T}}
) where T <: FieldElement
   R = parent(f)
   s = length(g)
   rcoeffs = T[]
   rexps = Vector{Int}[]
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

# check whether there exists a (possibly empty) w^'' such that
# w1 LM(g1) w2 = w1 LM(g1) w^'' LM(g2) u2
# and if that is the case, i.e. there is no overlap, returns false
# assumes that (w1, w2; u1, u2) are an obstruction of g1 and g2
# i.e. w1 LM(g1) w2 = u1 LM(g2) u2
function has_overlap(g1, g2, w1, w2, u1, u2)
    @assert vcat(w1, _leading_word(g1), w2) == vcat(u1, _leading_word(g2), u2)
    lw2 = _leading_word(g2)
    return length(w2) < length(lw2) + length(u2)
end

function is_redundant(
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

function get_obstructions(g::Vector{FreeAssAlgElem{T}}, s::Int) where T
   dummy_obstr = (Int[], Int[], Int[], Int[])
   empty_obstr_set = typeof(dummy_obstr)[]
   new_B = typeof(empty_obstr_set)[
              if i > j
                 empty_obstr_set
              elseif i == j
                 obstructions(_leading_word(g[i]))
              else
                 obstructions(_leading_word(g[i]), _leading_word(g[j]))
              end
           for i in 1:s, j in 1:s]
   obstr_count = 0
   for i in 1:s, j in 1:s
      obstr_count += length(new_B[i,j])
   end
   # TODO maybe here some redundancies can be removed too, check Kreuzer Xiu
   if groebner_debug_level > 0
      println("$obstr_count many obstructions")
   end
   return new_B
end


function add_obstructions(
   g::Vector{FreeAssAlgElem{T}},
   B::Matrix{Vector{NTuple{4, Vector{Int64}}}},
   s::Int
) where T
   dummy_obstr = (Int[], Int[], Int[], Int[])
   empty_obstr_set = typeof(dummy_obstr)[]
   new_B = Vector{NTuple{4, Vector{Int64}}}[
              if i < s && j < s
                 # copy old entry
                 B[i, j]
              elseif i > j
                 empty_obstr_set
              elseif i == j
                 obstructions(_leading_word(g[i]))
              else
                 obstructions(_leading_word(g[i]), _leading_word(g[j]))
              end
              for i in 1:s, j in 1:s]
   remove_redundancies!(new_B, s, g)
   return new_B
end


function groebner_basis_buchberger(
   g::Vector{FreeAssAlgElem{T}},
   degbound = typemax(Int)::Int
) where T <: FieldElement

   g = copy(g)
   R = parent(g[1])
   checked_obstructions = 0

   # step 1
   s = length(g)
   dummy_obstr = (Int[], Int[], Int[], Int[])
   empty_obstr_set = typeof(dummy_obstr)[]

   B = get_obstructions(g, s)
   while true
      @assert s == length(g)
      i = j = 0
      o = dummy_obstr
      # check all entries with i <= j
      for jj in 1:s, ii in 1:jj
         if !isempty(B[ii, jj])
           (i, j) = (ii, jj)
           o = popfirst!(B[i, j])
           break
         end
      end
      if !(i > 0 && j > 0)
         # B is empty
         return g
      end

      # step3
      S = _sub_rest(mul_term(inv(leading_coefficient(g[i])), o[1], g[i], o[2]),
                    mul_term(inv(leading_coefficient(g[j])), o[3], g[j], o[4]), 1)
      Sp = normal_form(S, g) # or normal_form_weak
      checked_obstructions += 1
      if groebner_debug_level > 0
         if checked_obstructions % 5000 == 0
            println("checked $checked_obstructions obstructions")
         end
      end
      if iszero(Sp) || total_degree(Sp) > degbound
         continue
      end

      # step4
      s += 1
      push!(g, Sp)
      if groebner_debug_level > 0
         println("adding new obstructions! checked $checked_obstructions so far")
      end
      B = add_obstructions(g, B, s)
   end
end

function groebner_basis(
   g::Vector{FreeAssAlgElem{T}},
   degbound = typemax(Int)::Int
) where T <: FieldElement
   return groebner_basis_buchberger(g, degbound)
end

