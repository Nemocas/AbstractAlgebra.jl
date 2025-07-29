pretty_lt(x, y) = false

pretty_lt(x::Integer, y::Integer) = isless(x, y)

pretty_lt(x::Rational{<:Integer}, y::Rational{<:Integer}) = isless(x, y)

# not sure about the default
pretty_eq(x, y) = (x == y)

function pretty_lt(x::T, y::T) where {T <: Union{Vector, Tuple}}
  return pretty_lt_it(x, y)
end

function pretty_eq(x::Vector, y::Vector)
  return all(pretty_eq.(x, y))
end

function pretty_lt_lex(x, y)
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
#  

struct PrettyOrderWrapper
  data::Any
end

Base.isless(x::PrettyOrderWrapper, y::PrettyOrderWrapper) = pretty_lt(x.data, y.data)

==(x::PrettyOrderWrapper, y::PrettyOrderWrapper) = pretty_eq(x.data, y.data)

pretty_sorting(x) = PrettyOrderWrapper(x)
