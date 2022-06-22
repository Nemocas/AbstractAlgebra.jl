#
#   FreeAssAhoCorasick.jl : implement bulk divide check for leading terms of free associative Algebra elements
#   for use e.g. in Groebner Basis computation
#
###############################################################################
#TODO how to properly export this?
export search, AhoCorasickAutomaton, AhoCorasickMatch, Word, insert_keyword!

using DataStructures

const Word = Vector{Int}

"""
Output stores for each node a tuple (i, k), where i is the index of the keyword k in 
the original list of keywords. If several keywords would be the output of the node, only
the one with the smallest index is stored
"""
mutable struct AhoCorasickAutomaton
    goto::Vector{Dict{Int, Int}}
    fail::Vector{Int}
    output::Vector{Tuple{Int, Word}}
end

struct AhoCorasickMatch
    last_position::Int
    keyword_index::Int
    keyword::Word
end

function AhoCorasickAutomaton(keywords::Vector{Word})
    automaton = AhoCorasickAutomaton([], [], [])
    construct_goto!(automaton, keywords)
    construct_fail!(automaton)
    return automaton
end

function lookup(automaton::AhoCorasickAutomaton, current_state::Int, next_letter::Int)
    ret_value = get(automaton.goto[current_state], next_letter, nothing)
    if current_state == 1 && isnothing(ret_value)
        return 1
    end
    return ret_value
    
end


function Base.length(automaton::AhoCorasickAutomaton)
    return length(automaton.goto)
end



function new_state!(automaton)
    push!(automaton.goto, Dict{Int, Int}())
    push!(automaton.output, (typemax(Int), []))
    push!(automaton.fail, 1)
    return length(automaton.goto)
end

function enter!(automaton::AhoCorasickAutomaton, keyword::Word, current_index)
    current_state = 1
        for c in keyword
        new_state = get(automaton.goto[current_state], c, nothing)
        if isnothing(new_state)
            new_state = new_state!(automaton)
            automaton.goto[current_state][c] = new_state
        end
        current_state = new_state
    end
    if automaton.output[current_state][1] > current_index
        automaton.output[current_state] = (current_index, keyword)
    end
end

function construct_goto!(automaton::AhoCorasickAutomaton, keywords::Vector{Word})
    new_state!(automaton)
    current_index = 1
    for keyword in keywords
        enter!(automaton, keyword, current_index)
        current_index += 1
    end
end

function construct_fail!(automaton::AhoCorasickAutomaton)
    q = Queue{Int}()
    for v in values(automaton.goto[1])
        enqueue!(q, v)
    end
    while !isempty(q)
        current_state = dequeue!(q)
        for k in keys(automaton.goto[current_state])
            new_state = lookup(automaton, current_state, k)
            enqueue!(q, new_state)
            state = automaton.fail[current_state]
            while isnothing(lookup(automaton, state, k))
                state = automaton.fail[state]
            end
            automaton.fail[new_state] = lookup(automaton, state, k)
            if automaton.output[new_state][1] > automaton.output[automaton.fail[new_state]][1]
                automaton.output[new_state] = automaton.output[automaton.fail[new_state]] # TODO check if this is the correct way to update output
            end

        end
    end
end

function insert_keyword!(aut::AhoCorasickAutomaton, keyword::Word, index::Int)
    enter!(aut, keyword, index)
    aut.fail = ones(Int, length(aut.goto))
    construct_fail!(aut)
end

"""

"""
function search(automaton::AhoCorasickAutomaton, word)
    current_state = 1
    result = AhoCorasickMatch(typemax(Int), typemax(Int), [])
    for i in 1:length(word)
        c = word[i]
        while true
            next_state = lookup(automaton, current_state, c)
            if !isnothing(next_state)
                current_state = next_state
                break
            else
                current_state = automaton.fail[current_state]
            end
        end
        if automaton.output[current_state][1] < result.keyword_index
            result = AhoCorasickMatch(i, automaton.output[current_state][1], automaton.output[current_state][2])
            #push!(output, (i, automaton.output[current_state]))
#            union!(output_set, (i, automaton.output[current_state]))
        end
    end
    if result.keyword == []
        return nothing
    end
    return result
end
