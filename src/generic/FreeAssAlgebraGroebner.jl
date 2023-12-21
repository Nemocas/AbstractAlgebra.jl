###############################################################################
#
#   FreeAssAlgebraGroebner.jl : free associative algebra R<x1,...,xn> Groebner basis
#
###############################################################################

include("PriorityQueue.jl")

const Monomial = Vector{Int}

abstract type Obstruction{T} end
"""
Represents the overlap of the leading term of the two polynomials
`first_poly` and `second_poly`. Here, `first_index` and `second_index`
are the indices of `first_poly` and `second_poly` respectively
in the Groebner basis.
The first and second prefix and suffix satisfy
first_prefix g_i first_suffix = second_prefix g_j second_suffix
where `g_i` is `first_poly` and `g_j` is `second_poly`.
"""
struct ObstructionTriple{T} <: Obstruction{T}
    first_poly::FreeAssAlgElem{T}
    second_poly::FreeAssAlgElem{T}
    first_prefix::Monomial
    first_suffix::Monomial
    second_prefix::Monomial
    second_suffix::Monomial
    first_index::Int
    second_index::Int
end

function ObstructionTriple{T}(first_poly::FreeAssAlgElem{T},
                              second_poly::FreeAssAlgElem{T},
                              pre_and_suffixes::NTuple{4, Monomial},
                              first_index::Int,
                              second_index::Int
                            ) where T
    return ObstructionTriple{T}(first_poly, second_poly, pre_and_suffixes[1], 
                                pre_and_suffixes[2], pre_and_suffixes[3], 
                                pre_and_suffixes[4], first_index, second_index)
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
    return FreeAssAlgElem{T}(parent(ot.first_poly), ot.first_prefix) *
           FreeAssAlgElem{T}(parent(ot.first_poly), _leading_word(ot.first_poly)) *
           FreeAssAlgElem{T}(parent(ot.first_poly), ot.first_suffix)
end

function s_polynomial(ot::ObstructionTriple{T}) where T
    first_term =
        FreeAssAlgElem{T}(parent(ot.first_poly), ot.first_prefix) *
        ot.first_poly *
        FreeAssAlgElem{T}(parent(ot.first_poly), ot.first_suffix)
    second_term =
        FreeAssAlgElem{T}(parent(ot.first_poly), ot.second_prefix) *
        ot.second_poly *
        FreeAssAlgElem{T}(parent(ot.first_poly), ot.second_suffix)
    return inv(leading_coefficient(ot.first_poly)) * first_term -
           inv(leading_coefficient(ot.second_poly)) * second_term
end

# skip all of the extra length-checking
function _leading_word(a::FreeAssAlgElem{T}) where T
    return a.exps[1]
end

@doc """
    gb_divides_leftmost(a::Word, aut::AhoCorasickAutomaton)

If an element of the Groebner basis that is stored in `aut` divides `a`,
return (true, a1, a2, keyword_index), where `keyword_index` is the index of the
keyword that divides `a` such that `a = a1 aut[keyword_index] a2`.
""" 
function gb_divides_leftmost(a::Word, aut::AhoCorasickAutomaton)
    match = search(aut, a)
    if isnothing(match)
        return (false, Word(), Word(), -1)
    end
    return (
        true,
        a[1:(match.last_position - length(match.keyword))],
        a[(match.last_position + 1):end],
        match.keyword_index,
    )
end

# implementation of the normal form function using aho corasick to check for all Groebner basis elements in parallel
@doc """
    normal_form(f::FreeAssAlgElem{T}, g::Vector{FreeAssAlgElem{T}}, aut::AhoCorasickAutomaton)

Assuming `g` is a Groebner basis and `aut` an Aho-Corasick automaton for the elements of `g`,
compute the normal form of `f` with respect to `g`
""" 
function normal_form(
    f::FreeAssAlgElem{T},
    g::Vector{FreeAssAlgElem{T}},
    aut::AhoCorasickAutomaton,
) where T
    R = parent(f)
    rexps = Monomial[]
    rcoeffs = T[]
    while length(f) > 0
        ok, left, right, match_index = gb_divides_leftmost(f.exps[1], aut)
        if ok
            qi = divexact(f.coeffs[1], g[match_index].coeffs[1])
            f = _sub_rest(f, mul_term(qi, left, g[match_index], right), 1)
        else
            push!(rcoeffs, f.coeffs[1])
            push!(rexps, f.exps[1])
            f = FreeAssAlgElem{T}(R, f.coeffs[2:end], f.exps[2:end], length(f) - 1)
        end
    end
    return FreeAssAlgElem{T}(R, rcoeffs, rexps, length(rcoeffs))
end

# normal form with leftmost word divisions
function normal_form(
    f::FreeAssAlgElem{T},
    g::Vector{FreeAssAlgElem{T}},
) where T<:FieldElement
    R = parent(f)
    s = length(g)
    rcoeffs = T[]
    rexps = Monomial[]
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
            f = FreeAssAlgElem{T}(R, f.coeffs[2:end], f.exps[2:end], length(f) - 1)
        end
    end
    r = FreeAssAlgElem{T}(R, rcoeffs, rexps, length(rcoeffs))
    return r
end

# weak normal form with leftmost word divisions
function normal_form_weak(
    f::FreeAssAlgElem{T},
    g::Vector{FreeAssAlgElem{T}},
) where T<:FieldElement
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

@doc raw"""
    interreduce!(g::Vector{FreeAssAlgElem{T}}) where T

Interreduce a given Groebner basis with itself, i.e. compute the normal form of each
element of `g` with respect to the rest of the elements and discard elements with
normal form $0$ and duplicates.
""" 
function interreduce!(g::Vector{FreeAssAlgElem{T}}) where T
    i = 1
    while length(g) > 1 && length(g) >= i
        aut = AhoCorasickAutomaton([g_j.exps[1] for g_j in g[1:end .!= i]])
        r = normal_form(g[i], g[1:end .!= i], aut)
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
function check_left_overlap(a::Monomial, b::Monomial, i::Int)
    if length(b) - i >= length(a)
        return false # this is a not a left overlap but might be a center overlap
    end
    for j in 0:(length(b) - i)
        if b[i + j] != a[j + 1]
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
function left_obstructions(a::Monomial, b::Monomial)
    v = Tuple{Monomial,Monomial}[]
    for i in 2:length(b)
        if check_left_overlap(a, b, i)
            if length(b) - i + 2 <= length(a) # w_2^' should not be empty!
                push!(v, (b[1:(i - 1)], a[(length(b) - i + 2):length(a)]))
            end
        end
    end
    return v
end


###
# find all non-trivial right-obstructions of a and b
# i.e. all words w_2 and w_1^' s.t. a w_1^' = w_2 b
# where length(w_1^') < length(b) and length(w_2) < length(a)
# the return vector is of the form [(w_1^', w_2), ...]
# if w_1^' or w_2 is empty, the corresponding obstruction is not returned
function right_obstructions(a::Monomial, b::Monomial)
    left_obstr = left_obstructions(b, a)
    right_obstr = [(w1, w2) for (w2, w1) in left_obstr]
    return right_obstr
end

###
# check whether a is a subword of b starting at index i
# a == b is also allowed
function check_center_overlap(a::Vector{Int}, b::Vector{Int}, i::Int)
    i + length(a) - 1 <= length(b) || return false
    return all(j -> a[j] == b[i + j - 1], 1:length(a))
end


function center_obstructions_first_in_second(a::Monomial, b::Monomial)
    v = Tuple{Monomial,Monomial}[]
    for i in 1:length(b)-length(a) + 1
        if check_center_overlap(a, b, i)
            push!(v, (b[1:(i - 1)], b[(i + length(a)):length(b)]))
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
function center_obstructions(a::Monomial, b::Monomial)
    if length(a) > length(b)
        return center_obstructions_first_in_second(b, a)
    else
        return center_obstructions_first_in_second(a, b)
    end
end

# all non-trivial ones
function obstructions(a::Monomial, b::Monomial)
    one = Int[] # the empty word
    res = Tuple{Monomial,Monomial,Monomial,Monomial}[]
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
function obstructions(a::Monomial)
    one = Int[] # the empty word
    res = Tuple{Monomial,Monomial,Monomial,Monomial}[]
    for x in left_obstructions(a, a)
        push!(res, (one, x[2], x[1], one))
    end
    for x in res
        @assert vcat(x[1], a, x[2]) == vcat(x[3], a, x[4])
    end
    return res
end


# check whether w_2 = v w_1 for some word v
function is_subword_right(w_1::Monomial, w_2::Monomial)
    if length(w_1) > length(w_2)
        return false
    end
    for i in 0:(length(w_1) - 1)
        if w_1[length(w_1) - i] != w_2[length(w_2) - i]
            return false
        end
    end
    return true
end


# check whether w_2 = w_1 v for some word v
function is_subword_left(w_1::Monomial, w_2::Monomial)
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
# check if obs2 is a subobstructon of obs1, i.e. if 
# the second pre-and suffix of obs2 are right- respectively left subwords of the second pre-and suffix of obs1.
# In other words, check if for obs1 = (w_i, w_i^'; u_j, u_j^') and obs2 = (w_k, w_k^'; v_l, v_l^')
# it holds that u_j == w v_l and u_j^' = v_l^' w^' for some w, w^'
# both w and w^' might be empty
function is_subobstruction(obs1_second_prefix::Monomial, obs1_second_suffix::Monomial, 
        obs2_second_prefix::Monomial, obs2_second_suffix)
    return is_subword_right(obs2_second_prefix, obs1_second_prefix) && is_subword_left(obs2_second_suffix, obs1_second_suffix)
end

function is_subobstruction(obs1::ObstructionTriple{T}, obs2::ObstructionTriple{T}) where T
    return is_subobstruction(obs1.second_prefix, obs1.second_suffix, obs2.second_prefix, obs2.second_suffix)
end

"""
if obs2 is a subobstruction of obs1, i.e. obs1[3] = w obs2[3] and obs1[4] = obs2[4]w',
returns length(ww')
thus, if it returns 0 and obs2 is a subobstruction of obs1, they are equal
if obs2 is not a subobstruction of obs1 the return value is useless
"""
function get_diff_length_for_subobstruction(
    obs1::ObstructionTriple{T},
    obs2::ObstructionTriple{T},
) where T
    return length(obs1.second_prefix) - length(obs2.second_prefix) +
           length(obs1.second_suffix) - length(obs2.second_suffix)
end

# check whether there exists a (possibly empty) w^'' such that
# w1 LM(g1) w2 = w1 LM(g1) w^'' LM(g2) u2
# only returns the correct value, if all the arguments come from an obstruction
# and if that is the case, i.e. there is no overlap, returns false
# assumes that (w1, w2; u1, u2) are an obstruction of g1 and g2
# i.e. w1 LM(g1) w2 = u1 LM(g2) u2
function has_overlap(g1, g2, w1, w2, u1, u2)
    @assert vcat(w1, _leading_word(g1), w2) == vcat(u1, _leading_word(g2), u2)
    lw2 = _leading_word(g2)
    return length(w2) < length(lw2) + length(u2)
end

function has_overlap(g2, w2, u2)
    lw2 = _leading_word(g2)
    return length(w2) < length(lw2) + length(u2)
end

function has_overlap(obs::ObstructionTriple{T}) where T
    return has_overlap(obs.second_poly, obs.first_suffix, obs.second_suffix)
end

function is_redundant(
    obs::ObstructionTriple{T},
    new_obstructions::PriorityQueue{Obstruction{T},FreeAssAlgElem{T}},
) where T
    # case 4b from Thm. 4.2.22 in Non-Commutative Groebner Bases and Applications by Xingqiang Xiu
    for obstruction_pair in new_obstructions
        o = obstruction_pair[1]
        if o.second_index == obs.second_index
            if is_subobstruction(obs, o)
                if get_diff_length_for_subobstruction(obs, o) > 0
                    return true
                elseif obs.first_index > o.first_index
                    return true
                elseif obs.first_index == o.first_index &&
                       get_diff_length_for_subobstruction(obs, o) == 0 &&
                       word_gt(obs.first_prefix, o.first_prefix)
                    return true
                end
            end
        end
    end
    return false
end

"""
check, whether obs1 is a proper multiple of obs2, i.e. they belong to the same polynomials and are of the form
obs1 = [w w_i, w_i' w'; w w_j, w_j' w'] and obs2 = [w_i, w_i'; w_j, w_j']
"""
function is_proper_multiple(
    obs1::ObstructionTriple{T},
    obs2::ObstructionTriple{T},
) where T
    if obs1.first_poly != obs2.first_poly || obs1.second_poly != obs2.second_poly #TODO compare indices instead?
        return false
    end
    if is_subword_right(obs2.first_prefix, obs1.first_prefix) &&
       is_subword_left(obs2.first_suffix, obs1.first_suffix)
        w = copy(obs1.first_prefix)
        w2 = copy(obs1.first_suffix)
        for _ in 1:length(obs2.first_prefix)
            pop!(w)
        end
        for _ in 1:length(obs2.first_suffix)
            popfirst!(w2)
        end
        if length(w) + length(w2) == 0
            return false
        end
        @assert obs1.first_prefix == vcat(w, obs2.first_prefix)
        @assert obs1.first_suffix == vcat(obs2.first_suffix, w2)
        return obs1.second_prefix == vcat(w, obs2.second_prefix) &&
               obs1.second_suffix == vcat(obs2.second_suffix, w2)
    else
        return false
    end
end

"""
check, whether obs is a proper multiple of any of the obstructions in the priority queue
"""
function is_proper_multiple(
    obs::ObstructionTriple{T},
    obstructions::PriorityQueue{Obstruction{T},FreeAssAlgElem{T}},
) where T
    for obspair in obstructions
        obs2 = obspair[1]
        if is_proper_multiple(obs, obs2)
            return true
        end
    end
    return false
end

function is_redundant(
    obs::ObstructionTriple{T},
    new_obstructions::PriorityQueue{Obstruction{T},FreeAssAlgElem{T}},
    newest_element::FreeAssAlgElem{T},
    newest_index::Int,
) where T
    w1 = []
    w2 = []
    for i in 1:length(obs.second_poly.exps[1])
        word_to_check =
            vcat(obs.second_prefix, obs.second_poly.exps[1], obs.second_suffix)
        if check_center_overlap(newest_element.exps[1], word_to_check, i)
            w1 = word_to_check[1:(i - 1)]
            w2 = word_to_check[(i + length(newest_element.exps[1])):end]
            break
        end
    end
    if length(w1) + length(w2) == 0
        return false
    end
    obs1 = ObstructionTriple{T}(
        obs.first_poly,
        newest_element,
        obs.first_prefix, obs.first_suffix, w1, w2,
        obs.first_index,
        newest_index,
    )
    obs2 = ObstructionTriple{T}(
        obs.second_poly,
        newest_element,
        obs.second_prefix, obs.second_suffix, w1, w2,
        obs.second_index,
        newest_index,
    )
    o1_bool = !has_overlap(obs1) || is_proper_multiple(obs1, new_obstructions) # TODO maybe only call is_proper_multiple if both obs have no overlap for performance?
    o2_bool = !has_overlap(obs2) || is_proper_multiple(obs2, new_obstructions)
    return o1_bool && o2_bool
end


function remove_redundancies!(
    all_obstructions::PriorityQueue{Obstruction{T},FreeAssAlgElem{T}},
    newest_index::Int,
    newest_element::FreeAssAlgElem{T},
) where T
    del_counter = 0
    new_obstructions = PriorityQueue{Obstruction{T},FreeAssAlgElem{T}}()
    old_obstructions = PriorityQueue{Obstruction{T},FreeAssAlgElem{T}}()
    for obstr_pair in all_obstructions
        if obstr_pair[1].second_index == newest_index
            new_obstructions[obstr_pair[1]] = obstr_pair[2]
        else
            old_obstructions[obstr_pair[1]] = obstr_pair[2]
        end
    end
    for obstr_pair in new_obstructions
        if is_redundant(obstr_pair[1], new_obstructions)
            del_counter += 1
            delete!(new_obstructions, obstr_pair[1])
            delete!(all_obstructions, obstr_pair[1])
        end
    end
    if length(new_obstructions) == 0
        return nothing
    end

    for obstr_pair in old_obstructions
        if is_redundant(obstr_pair[1], new_obstructions, newest_element, newest_index)
            del_counter += 1
            delete!(all_obstructions, obstr_pair[1])
        end
    end
    # TODO case 4e from Thm 4.1 in Kreuzer Xiu
end

function get_obstructions(g::Vector{FreeAssAlgElem{T}}) where T
    s = length(g)
    result = PriorityQueue{Obstruction{T},FreeAssAlgElem{T}}()
    for i in 1:s, j in 1:i
        if i == j
            obs = obstructions(_leading_word(g[i]))
        else
            obs = obstructions(_leading_word(g[i]), _leading_word(g[j]))
        end
        for o in obs
            triple = ObstructionTriple{T}(g[i], g[j], o, i, j)
            push!(result, triple => common_multiple_leading_term(triple))
        end
    end
    # TODO maybe here some redundancies can be removed too, check Kreuzer Xiu
    return result
end


function add_obstructions!(
    obstruction_queue::PriorityQueue{Obstruction{T},FreeAssAlgElem{T}},
    g::Vector{FreeAssAlgElem{T}},
) where T
    s = length(g)
    for i in 1:s
        if i == s
            obs = obstructions(_leading_word(g[i]))
        else
            obs = obstructions(_leading_word(g[i]), _leading_word(g[s]))
        end
        for o in obs
            triple = ObstructionTriple{T}(g[i], g[s], o, i, s)
            push!(obstruction_queue, triple => common_multiple_leading_term(triple))
        end
    end
    #remove_redundancies!(obstruction_queue, s, g[s]) #TODO too slow in practice
end


function groebner_basis_buchberger(
    g::Vector{FreeAssAlgElem{T}},
    reduction_bound::Int = typemax(Int),
    remove_redundancies::Bool = false
) where T<:FieldElement
    g = copy(g)
    #   interreduce!(g) # on some small examples, this increases running time, so it might not be optimal to use this here
    nonzero_reductions = 0
    # compute the aho corasick automaton
    # to make normal form computation more efficient
    aut = AhoCorasickAutomaton([g_i.exps[1] for g_i in g])
    # step 1 from Thm. 5.2.12 Noncommutative Groebner Bases and Applications, Xingqiang Xiu
    obstruction_queue = get_obstructions(g)
    while !isempty(obstruction_queue) # step 2
        obstruction = popfirst!(obstruction_queue)[1]
        # step3
        S = s_polynomial(obstruction)
        Sp = normal_form(S, g, aut) # or normal_form_weak
        if iszero(Sp)
            continue
        end
        nonzero_reductions += 1
        # step 4
        push!(g, Sp)
        insert_keyword!(aut, Sp.exps[1], length(g))
        if nonzero_reductions >= reduction_bound
            return g
        end
        add_obstructions!(obstruction_queue, g)
        if remove_redundancies
            remove_redundancies!(obstruction_queue, length(g), g[length(g)])
        end
    end
    return g
end

@doc """
    groebner_basis(g::Vector{FreeAssAlgElem{T}}, reduction_bound::Int = typemax(Int), remove_redundancies::Bool = false)

Compute a Groebner basis for the ideal generated by `g`. Stop when `reduction_bound` many
non-zero entries have been added to the Groebner basis. If the computation stops due to the bound being exceeded, 
the result is in general not an actual Groebner basis, just a subset of one. However, whenever the normal form with
respect to this incomplete Groebner basis is `0`, it will also be `0` with respect to the full Groebner basis.
""" 
function groebner_basis(
    g::Vector{FreeAssAlgElem{T}},
    reduction_bound::Int = typemax(Int),
    remove_redundancies::Bool = false
) where T<:FieldElement
    return groebner_basis_buchberger(g, reduction_bound, remove_redundancies)
end
