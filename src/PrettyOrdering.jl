module PrettyOrdering

pretty_lt(x, y) = false

pretty_lt(x::Integer, y::Integer) = isless(x, y)

pretty_lt(x::Rational{<:Integer}, y::Rational{<:Integer}) = isless(x, y)

# not sure about the default
pretty_eq(x, y) = (x === y)

function pretty_lt(x::T, y::S) where {T <: Union{Vector, Tuple}, S <: Union{Vector, Tuple}}
  return pretty_lt_lex(x, y)
end

function pretty_eq(x::T, y::S) where {T <: Union{Vector, Tuple}, S <: Union{Vector, Tuple}}
  length(x) == length(y) || return false
  return all(pretty_eq.(x, y))
end

function pretty_lt_lex(x, y)
  length(x) < length(y) && return true
  length(x) > length(y) && return false
  # assume both things have same length
  for (a, b) in zip(x, y)
    if pretty_eq(a, b)
      continue
    end
    if pretty_lt(a, b)
      return true
    else
      return false
    end
  end
  # pretty_eq
  return false
end

################################################################################
#
#  Wrapper type
#
################################################################################

struct PrettyOrderWrapper{T}
  data::T

  PrettyOrderWrapper(x) = new{typeof(x)}(x)
end

Base.isless(x::PrettyOrderWrapper, y::PrettyOrderWrapper) = pretty_lt(x.data, y.data)

Base.:(==)(x::PrettyOrderWrapper, y::PrettyOrderWrapper) = pretty_eq(x.data, y.data)

pretty_sort(x) = PrettyOrderWrapper(x)

end # module

import ..PrettyOrdering: pretty_eq
import ..PrettyOrdering: pretty_lt
import ..PrettyOrdering: pretty_lt_lex
import ..PrettyOrdering: pretty_sort
