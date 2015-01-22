type PariVector{T <: PariRing}
   base_ring::PariRing
end

type pari_vec{T <: PariRing}
   data::Ptr{Int}
   parent::PariVector{T}

   function pari_vec{R <: PariRing}(v::Ptr{Int}, par::R)
      r = new(gclone(v), PariVector{T}(par))
      finalizer(r, _pari_vec_unclone)
      return r
   end
end

_pari_vec_unclone(a::pari_vec) = gunclone(a.data)

function getindex(a::pari_vec, n::Int)
   a.parent.base_ring(reinterpret(Ptr{Int}, unsafe_load(a.data + n*sizeof(Int))))
end