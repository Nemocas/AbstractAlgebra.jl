###############################################################################
#
#   FreeAssAhoCorasick.jl : implement bulk divide check for leading terms of free associative Algebra elements
#   for use e.g. in Groebner Basis computation
#
###############################################################################
using DataStructures

const Word = Vector{Int}

struct AhoCorasickAutomaton
    goto::Vector{Dict{Int, Int}}
    fail::Vector{Int}
    output::Vector{Set{Word}}
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
    push!(automaton.output, Set{Word}())
    push!(automaton.fail, 1)
    return length(automaton.goto)
end

function enter!(automaton::AhoCorasickAutomaton, keyword::Word)
    current_state = 1
    for c in keyword
        new_state = get(automaton.goto[current_state], c, nothing)
        if isnothing(new_state)
            new_state = new_state!(automaton)
            automaton.goto[current_state][c] = new_state
        end
        current_state = new_state
    end

    push!(automaton.output[current_state], keyword)

end

function construct_goto!(automaton::AhoCorasickAutomaton, keywords::Vector{Word})
    new_state!(automaton)
    for keyword in keywords
        enter!(automaton, keyword)
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
            union!(automaton.output[new_state], automaton.output[automaton.fail[new_state]])
        end
    end
end

"""

"""
function search(automaton::AhoCorasickAutomaton, word)
    current_state = 1
    output_set = Set{Word}()
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
        if !isempty(automaton.output[current_state])
            println(i)
        end
        union!(output_set, automaton.output[current_state])
    end
    return output_set
end
