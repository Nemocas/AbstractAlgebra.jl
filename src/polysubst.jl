for T in subtypes(PolyRingElem)
  (f::T)(a) = subst(f, a)

  function (f::T)(a::T)
     if parent(f) != parent(a)
        return subst(f, a)
     end
     return compose(f, a; inner = :second)
  end

  (f::T)(a::Integer) = evaluate(f, a)

  function (f::T)(a::RingElem)
     if parent(a) != base_ring(f)
        return subst(f, a)
     end
     return evaluate(f, a)
  end
end

for T in subtypes(NCPolyRingElem)
  (f::T)(a::Integer) = evaluate(f, a)

  function (f::T)(a::NCRingElem)
     return evaluate(f, a)
  end
end
