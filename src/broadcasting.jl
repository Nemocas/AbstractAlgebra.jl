################################################################################
#
#  Broadcasting for MatrixElem
#
################################################################################

################################################################################
#
#  Broadcaststyle
#
################################################################################

# We would like to use Broadcast.ArrayStyle{MatElem}, as our "broadcasting
# style" (whatever this means), but this is not allowed since
# !(MatElem <: AbstractArray). We work around this by using a dummy type.
# This will only be used in the similar method, where we can just ignore it.

struct BroadcastDummy <: AbstractMatrix{Int}
end

Base.broadcastable(x::MatElem) = x

Base.BroadcastStyle(::Type{<:MatElem}) = Broadcast.ArrayStyle{BroadcastDummy}()

################################################################################
#
#  Creation of destination matrix
#
################################################################################

# I use a non-recursive version at the leaft to not confusing inference too much
function _compute_target_size(a::Tuple)
  if length(a) == 1
    return _compute_target_size_nonrec(a[1])
  else
    s = _compute_target_size_nonrec(a[1])
    return _promote_size(s, _compute_target_size(Base.tail(a)))
  end
end

function _promote_size(s1, s2)
  if s1 isa Tuple{}
    return s2
  end
  return s1
end

function _compute_target_size_nonrec(a::Broadcast.Broadcasted)
  return _compute_target_size(a.args)
end

function _compute_target_size_nonrec(a::MatElem)
  return size(a)
end

function _compute_target_size_nonrec(a::NCRingElement)
  return ()
end

function _compute_target_size_nonrec(a::Ref)
  return ()
end

function _compute_target_size_nonrec(a)
  error("Broadcasting for type $(typeof(a)) not implemented")
end

function _promote_dest_func(f::T, a::Tuple) where {T}
  g = _compute_elements(a)
  K = parent(f(g...))
  s = _compute_target_size(a)
  return zero_matrix(K, s...)
end

# Same as above about the leafs
function _compute_elements(a::Tuple)
  if length(a) == 1
    return (_compute_elements_nonrec(a[1]),)
  else
    return (_compute_elements_nonrec(a[1]), _compute_elements(Base.tail(a))...)
  end
end

function _compute_elements_nonrec(a::Broadcast.Broadcasted)
  return _promote_dest_func_elem(a.f, a.args)
end

function _compute_elements_nonrec(a::MatElem)
  return zero(base_ring(a))
end

function _compute_elements_nonrec(a::NCRingElement)
  return a
end

function _compute_elements_nonrec(a::Ref)
  return a.x
end

function _compute_elements_nonrec(a)
  error("Broadcasting for type $(typeof(a)) not implemented")
end

function _promote_dest_func_elem(f::T, a::Tuple) where {T}
  g = _compute_elements(a)
  K = parent(f(g...))
  return zero(K)
end

################################################################################
#
#  Similar functionality
#
################################################################################

# The most tricky part:
#
# Any expression of the form fun.(args...) is transformed into a
# Broadcasted object bc with bc.fun == and bc.args == x (a Tuple)
# Note that the x can themselves be again objects of type Broadcasted.
# So we end up with a nice expression tree and our aim is to find
# the shape and coefficient ring of the final output matrix
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{BroadcastDummy}}, ::Type{ElType}) where ElType
  dest = _promote_dest_func(bc.f, bc.args)
  return similar(dest, (length.(axes(bc)))...)
end

Base.copyto!(dest::MatElem, bc::Broadcast.Broadcasted) = Base.copyto!(dest, convert(Broadcast.Broadcasted{Nothing}, bc))

function Base.copyto!(dest::MatElem, bc::Broadcast.Broadcasted{Nothing})
  axes(dest) == axes(bc) || throwdm(axes(dest), axes(bc))
  # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
  if bc.f === identity && bc.args isa Tuple{AbstractArray} # only a single input argument to broadcast!
    A = bc.args[1]
    if axes(dest) == axes(A)
      return copyto!(dest, A)
    end
  end
  bcc = Broadcast.preprocess(dest, bc)
  @inbounds for I in eachindex(bcc)
    dest[I] = bcc[I]
  end
  return dest
end
