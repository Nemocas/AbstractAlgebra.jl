import Base.hash
export hash

################################
# fun with fmpz ################
#
# should of course be done in c
################################

function fmpz_is_small(a::fmpz)
  return _fmpz_is_small(a.d)
end

function fmpz_limbs(a::fmpz)
  return fmpz_limbs(a.d)
end

function hash_integer(a::fmpz, h::UInt)  ## cloned from the BigInt function
                                         ## in Base
  return _hash_integer(a.d, h)
end

function hash(a::fmpz, h::UInt)
#  pow = trailing_zeros(a)
#  h = Base.hash_uint(unsigned(pow) $ h) $ h
#  h = hash_integer(a>>pow, h)
  h = hash_integer(a, h)
  return h
end

function _fmpz_is_small(a::Int)
  return (unsigned(a)>>(WORD_SIZE-2) !=1)
end

function fmpz_limbs(a::Int)
  if _fmpz_is_small(a)
    return 0
  end
  b = unsafe_load(convert(Ptr{Cint}, unsigned(a)<<2), 2)
  return b
end

function _hash_integer(a::Int, h::UInt)  ## cloned from the BigInt function
                                         ## in Base
  s = fmpz_limbs(a)
  s == 0 && return Base.hash_integer(a, h)
  # get the pointer after the first two Cint
  d = convert(Ptr{Ptr{UInt}}, unsigned(a)<<2) + 2sizeof(Cint)
  p = unsafe_load(d)
  b = unsafe_load(p)
  h = Base.hash_uint(ifelse(s < 0, -b, b) $ h) $ h
  for k = 2:abs(s)
    h = Base.hash_uint(unsafe_load(p, k) $ h) $ h
  end
  return h
end

function hash(a::fmpq, h::UInt)
  return _hash_integer(a.num, _hash_integer(a.den, h))
end

