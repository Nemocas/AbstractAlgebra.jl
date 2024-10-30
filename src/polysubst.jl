(f::PolyRingElem)(a) = subst(f, a)

function (f::PolyRingElem)(a::PolyRingElem)
    typeof(f) == typeof(a) || return subst(f, a)
    parent(f) == parent(a) || return subst(f, a)
    return compose(f, a; inner = :second)
end

(f::PolyRingElem)(a::Integer) = evaluate(f, a)

function (f::PolyRingElem)(a::RingElem)
    base_ring(f) == parent(a) || return subst(f, a)
    return evaluate(f, a)
end

(f::NCPolyRingElem)(a::Integer) = evaluate(f, a)

function (f::NCPolyRingElem)(a::NCRingElem)
    return evaluate(f, a)
end
