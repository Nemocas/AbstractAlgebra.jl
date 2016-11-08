###############################################################################
#
#   SRPoly.jl : Generic recursive sparse multivariate polynomials over rings
#
###############################################################################

export GenSRPoly, GenSRPolyRing, SRPolynomialRing

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type{T, K}(::Type{GenSRPoly{T, K}}) = GenSRPolyRing{T, K}

elem_type{T <: RingElem, K <: RingElem}(::GenSRPolyRing{T, K}) = GenSRPoly{T, K}

var(a::GenSRPolyRing) = a.S

function gen{T, K <:RingElem}(a::GenSRPolyRing{T, K})
   return a([one(base_ring(a))], [UInt(1)])
end

function gens{T <: RingElem, K <: RingElem}(a::GenSRPolyRing{T, K})
   A = Array(GenSRPoly, a.num_vars)
   R = a
   for i = 1:a.num_vars
      A[i] = gen(R)
      R = base_ring(R)
   end
   return A
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function coeff(x::GenSRPoly, i::Int)
   i < 0 && throw(DomainError())
   return x.coeffs[i + 1]
end

num_vars(x::GenSRPoly) = parent(x).num_vars

one(R::GenSRPolyRing) = R(1)

zero(R::GenSRPolyRing) = R(0)
 
iszero(x::GenSRPoly) = length(x) = 0

isone(x::GenSRPoly) = x == 1

length(x::GenSRPoly) = x.length

lead(x::GenSRPoly) = x.length == 0 ? base_ring(x)() : x.coeffs[x.length]

function normalise(a::GenSRPoly, n::Int)
   while n > 0 && iszero(a.coeffs[n]) 
      n -= 1
   end
   return n
end

base_ring(a::GenSRPoly) = base_ring(parent(a))

base_ring{T <: RingElem, K <: RingElem}(R::GenSRPolyRing{T, K}) = R.base_ring::parent_type(T)

parent(a::GenSRPoly) = a.parent

###############################################################################
#
#   String I/O
#
###############################################################################

function show{T <: RingElem, K <: RingElem}(io::IO, x::GenSRPoly{T, K})
    len = length(x)
    U = string(var(parent(x)))
    if len == 0
      print(io, base_ring(x)(0))
    else
      for i = 1:len
        c = coeff(x, len - i)
        bracket = needs_parentheses(c)
        if i != 1 && !is_negative(c)
          print(io, "+")
        end
        X = x.exps[len - i + 1]
        if !isone(c) && (c != -1 || show_minus_one(typeof(c)))
          if bracket
            print(io, "(")
          end
          show(io, c)
          if bracket
            print(io, ")")
          end
          if c != 1 && !(c == -1 && !show_minus_one(typeof(c))) && X != UInt(0)
             print(io, "*")
          end
        end
        if c == -1 && !show_minus_one(typeof(c))
          print(io, "-")
        end
        if X == UInt(0)
          if c == 1
             print(io, c)
          elseif c == -1 && !show_minus_one(typeof(c))
             print(io, 1)
          end
        end
        n = reinterpret(Int, X)
        if n != 0
           print(io, U)
        end
        if n > 1
           print(io, "^", n)
        end
      end
    end
end

function show(io::IO, p::GenSRPolyRing)
   const max_vars = 5 # largest number of variables to print
   n = p.num_vars
   print(io, "Recursive Multivariate Polynomial Ring in ")
   if n > max_vars
      print(io, p.num_vars)
      print(io, " variables ")
   end
   for i = 1:min(n - 1, max_vars - 1)
      print(io, string(p.S), ", ")
      p = base_ring(p)
   end
   if n > max_vars
      print(io, "..., ")
   end
   print(io, string(p.S))
   print(io, " over ")
   show(io, base_ring(p))
end

show_minus_one{T <: RingElem, K <: RingElem}(::Type{GenSRPoly{T, K}}) = show_minus_one(T)

needs_parentheses{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}) = length(a) > 1

is_negative(x::GenSRPoly) = length(x) <= 1 && is_negative(coeff(x, 0))

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function -{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K})
   r = parent(a)()
   fit!(r, length(a))
   for i = 1:length(a)
      r.coeffs[i] = -a.coeffs[i]
      r.exps[i] = a.exps[i]
   end
   r.length = a.length
   return r
end

function +{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   par = parent(a)
   r = par()
   fit!(r, length(a) + length(b))
   i = 1
   j = 1
   k = 1
   while i <= length(a) && j <= length(b)
      if a.exps[i] < b.exps[j]
         r.coeffs[k] = a.coeffs[i]
         r.exps[k] = a.exps[i]
         i += 1
      elseif a.exps[i] == b.exps[j]
         c = a.coeffs[i] + b.coeffs[j]
         if c != 0
            r.coeffs[k] = c
            r.exps[k] = a.exps[i]
         else
            k -= 1
         end
         i += 1
         j += 1
      else
         r.coeffs[k] = b.coeffs[j]
         r.exps[k] = b.exps[j]
         j += 1
      end
      k += 1
   end
   while i <= length(a)
      r.coeffs[k] = a.coeffs[i]
      r.exps[k] = a.exps[i]
      i += 1
      k += 1
   end
   while j <= length(b)
      r.coeffs[k] = b.coeffs[j]
      r.exps[k] = b.exps[j]
      j += 1
      k += 1
   end
   r.length = k - 1
   return r
end

function -{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   par = parent(a)
   r = par()
   fit!(r, length(a) + length(b))
   i = 1
   j = 1
   k = 1
   while i <= length(a) && j <= length(b)
      if a.exps[i] < b.exps[j]
         r.coeffs[k] = a.coeffs[i]
         r.exps[k] = a.exps[i]
         i += 1
      elseif a.exps[i] == b.exps[j]
         c = a.coeffs[i] - b.coeffs[j]
         if c != 0
            r.coeffs[k] = c
            r.exps[k] = a.exps[i]
         else
            k -= 1
         end
         i += 1
         j += 1
      else
         r.coeffs[k] = -b.coeffs[j]
         r.exps[k] = b.exps[j]
         j += 1
      end
      k += 1
   end
   while i <= length(a)
      r.coeffs[k] = a.coeffs[i]
      r.exps[k] = a.exps[i]
      i += 1
      k += 1
   end
   while j <= length(b)
      r.coeffs[k] = -b.coeffs[j]
      r.exps[k] = b.exps[j]
      j += 1
      k += 1
   end
   r.length = k - 1
   return r
end

immutable heap_sr
   exp::UInt
   n::Int
end

function isless(a::heap_sr, b::heap_sr)
   return a.exp < b.exp
end

function ==(a::heap_sr, b::heap_sr)
   return a.exp == b.exp
end

function heapinsert!(xs::Array{heap_sr, 1}, ys::Array{heap_t, 1}, m::Int, exp::UInt)
   i = n = length(xs) + 1
   @inbounds if i != 1 && exp == xs[1].exp
      ys[m] = heap_t(ys[m].i, ys[m].j, xs[1].n)
      xs[1] = heap_sr(exp, m)
      return
   end
   @inbounds while (j = heapparent(i)) >= 1
      if exp == xs[j].exp
         ys[m] = heap_t(ys[m].i, ys[m].j, xs[j].n)
         xs[j] = heap_sr(exp, m)
         return
      elseif exp < xs[j].exp
         i = j
      else
         break
      end
   end
   push!(xs, heap_sr(exp, 0))
   @inbounds while n > i
      xs[n] = xs[heapparent(n)]
      n >>= 1
   end
   xs[i] = heap_sr(exp, m)
   return
end

function heappop!(xs::Array{heap_sr, 1})
   s = length(xs)
   x = xs[1]
   i = 1
   j = 2
   @inbounds while j < s
      if xs[j].exp >= xs[j + 1].exp
         j += 1
      end
      xs[i] = xs[j]
      i = j
      j *= 2
   end
   exp = xs[s].exp
   j = i >> 1
   @inbounds while i > 1 && exp < xs[j].exp
      xs[i] = xs[j]
      i = j
      j >>= 1
   end
   xs[i] = xs[s]
   pop!(xs)
   return
end

function mul_johnson{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   par = parent(a)
   R = base_ring(par)
   m = length(a)
   n = length(b)
   if m == 0 || n == 0
      return par()
   end
   H = Array(heap_sr, 0)
   I = Array(heap_t, 0)
   # set up heap
   push!(H, heap_sr(a.exps[1] + b.exps[1], 1))
   push!(I, heap_t(1, 1, 0))
   r_alloc = max(m, n) + n
   Rc = Array(T, r_alloc)
   Re = Array(UInt, r_alloc)
   k = 0
   c = R()
   Q = Array(Int, 0)
   while !isempty(H)
      exp = H[1].exp
      k += 1
      if k > r_alloc
         r_alloc *= 2
         resize!(Rc, r_alloc)
         resize!(Re, r_alloc)
      end
      first = true
      @inbounds while !isempty(H) && H[1].exp == exp
         x = H[1]
         heappop!(H)
         v = I[x.n]
         if first
            Rc[k] = a.coeffs[v.i]*b.coeffs[v.j]
            Re[k] = exp
            first = false
         else
            addmul!(Rc[k], a.coeffs[v.i], b.coeffs[v.j], c)
         end
         if v.j < n || v.j == 1
            push!(Q, x.n)
         end
         while (xn = v.next) != 0
            v = I[xn]
            addmul!(Rc[k], a.coeffs[v.i], b.coeffs[v.j], c)
            if v.j < n || v.j == 1
               push!(Q, xn)
            end
         end
      end
      @inbounds while !isempty(Q)
         xn = pop!(Q)
         v = I[xn]
         if v.j == 1 && v.i < m
            push!(I, heap_t(v.i + 1, 1, 0))
            Collections.heappush!(H, heap_sr(a.exps[v.i + 1] + b.exps[1], length(I)))       
         end
         if v.j < n
            I[xn] = heap_t(v.i, v.j + 1, 0)
            heapinsert!(H, I, xn, a.exps[v.i] + b.exps[v.j + 1]) # either chain or insert v into heap   
         end
      end
      if Rc[k] == 0
         k -= 1
      end
   end
   resize!(Rc, k)
   resize!(Re, k)
   return parent(a)(Rc, Re)
end

function *{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   m = length(a)
   n = length(b)
   if m == 0 || n == 0
      return parent(a)()
   end
   if m < n
      return b*a
   end
   r = a*b.coeffs[1]
   shift_left!(r, b.exps[1])
   for i = 2:n
      s = a*b.coeffs[i]
      shift_left!(s, b.exps[i])
      r += s
   end
   return r
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   par = parent(a)
   R = base_ring(par)
   m = length(a)
   n = length(b)
   maxn = a.exps[m]
   n == 0 && error("Division by zero in divrem_monagan_pearce")
   if a.exps[m] < b.exps[n]
      return par(), a
   end
   H = Array(heap_sr, 0)
   I = Array(heap_t, 0)
   # set up heap
   push!(H, heap_sr(maxn - a.exps[m], 1))
   push!(I, heap_t(0, 1, 0))
   q_alloc = max(m - n, n)
   r_alloc = n
   Qc = Array(T, q_alloc)
   Qe = Array(UInt, q_alloc)
   Rc = Array(T, r_alloc)
   Re = Array(UInt, r_alloc)
   k = 0
   l = 0
   s = n
   c = R()
   qc = R()
   m1 = -R(1)
   mb = -b.coeffs[n]
   Q = Array(Int, 0)
   reuse = Array(Int, 0)
   while !isempty(H)
      exp = H[1].exp
      k += 1
      if k > q_alloc
         q_alloc *= 2
         resize!(Qc, q_alloc)
         resize!(Qe, q_alloc)
      end
      @inbounds while !isempty(H) && H[1].exp == exp
         x = H[1]
         heappop!(H)
         v = I[x.n]
         if v.i == 0
            addmul!(qc, a.coeffs[m + 1 - v.j], m1, c)
         else
            addmul!(qc, b.coeffs[n + 1 - v.i], Qc[v.j], c)
         end
         if v.i != 0 || v.j < m
            push!(Q, x.n)
         else
            push!(reuse, x.n)
         end
         while (xn = v.next) != 0
            v = I[xn]
            if v.i == 0
               addmul!(qc, a.coeffs[m + 1 - v.j], m1, c)
            else
               addmul!(qc, b.coeffs[n + 1 - v.i], Qc[v.j], c)
            end
            if v.i != 0 || v.j < m
               push!(Q, xn)
            else
               push!(reuse, xn)
            end
         end
      end
      @inbounds while !isempty(Q)
         xn = pop!(Q)
         v = I[xn]
         if v.i == 0
            I[xn] = heap_t(0, v.j + 1, 0)
            heapinsert!(H, I, xn, maxn - a.exps[m - v.j]) # either chain or insert into heap   
         elseif v.j < k - 1
            I[xn] = heap_t(v.i, v.j + 1, 0)
            heapinsert!(H, I, xn, maxn - b.exps[n + 1 - v.i] - Qe[v.j + 1]) # either chain or insert into heap
         elseif v.j == k - 1
            s += 1
            push!(reuse, xn)
         end  
      end
      if qc == 0
         k -= 1
      else
         d1 = maxn - b.exps[n] >= exp
         texp = maxn - b.exps[n] - exp
         if !d1
            l += 1
            if l >= r_alloc
               r_alloc *= 2
               resize!(Rc, r_alloc)
               resize!(Re, r_alloc)
            end
            Rc[l] = -qc
            Re[l] = maxn - exp
            k -= 1
         else
            tq, tr = divrem(qc, mb)
            if tr != 0
               l += 1
               if l >= r_alloc
                  r_alloc *= 2
                  resize!(Rc, r_alloc)
                  resize!(Re, r_alloc)
               end
               Rc[l] = -tr
               Re[l] = maxn - exp 
            end
            if tq != 0
               Qc[k] = tq
               Qe[k] = texp
               for i = 2:s
                  if !isempty(reuse)
                     xn = pop!(reuse)
                     I[xn] = heap_t(i, k, 0)
                     heapinsert!(H, I, xn, maxn - b.exps[n + 1 - i] - Qe[k]) # either chain or insert into heap
                  else
                     push!(I, heap_t(i, k, 0))
                     Collections.heappush!(H, heap_sr(maxn - b.exps[n + 1 - i] - Qe[k], length(I)))
                  end
               end                 
               s = 1
            else
               k -= 1
            end
         end
      end
      zero!(qc)
   end
   resize!(Qc, k)
   resize!(Qe, k)
   resize!(Rc, l)
   resize!(Re, l)
   reverse!(Qc)
   reverse!(Qe)
   reverse!(Rc)
   reverse!(Re)
   return parent(a)(Qc, Qe), parent(a)(Rc, Re)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

function *{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, n::Integer)
   r = parent(a)()
   fit!(r, length(a))
   j = 1
   for i = 1:length(a)
      c = a.coeffs[i]*n
      if c != 0
         r.coeffs[j] = c 
         r.exps[j] = a.exps[i]
         j += 1
      end
   end
   r.length = j - 1
   return r
end

function *{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, n::fmpz)
   r = parent(a)()
   fit!(r, length(a))
   j = 1
   for i = 1:length(a)
      c = a.coeffs[i]*n
      if c != 0
         r.coeffs[j] = c 
         r.exps[j] = a.exps[i]
         j += 1
      end
   end
   r.length = j - 1
   return r
end

function *{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, n::T)
   r = parent(a)()
   fit!(r, length(a))
   j = 1
   for i = 1:length(a)
      c = a.coeffs[i]*n
      if c != 0
         r.coeffs[j] = c 
         r.exps[j] = a.exps[i]
         j += 1
      end
   end
   r.length = j - 1
   return r
end

*{T <: RingElem, K <: RingElem}(n::T, a::GenSRPoly{T, K}) = a*n

*{T <: RingElem, K <: RingElem}(n::Integer, a::GenSRPoly{T, K}) = a*n

*{T <: RingElem, K <: RingElem}(n::fmpz, a::GenSRPoly{T, K}) = a*n

+{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::Integer) = a + parent(a)(b)

+{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::fmpz) = a + parent(a)(b)

-{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::Integer) = a - parent(a)(b)

-{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::fmpz) = a - parent(a)(b)

+{T <: RingElem, K <: RingElem}(a::Integer, b::GenSRPoly{T, K}) = parent(b)(a) + b

+{T <: RingElem, K <: RingElem}(a::fmpz, b::GenSRPoly{T, K}) = parent(b)(a) + b

-{T <: RingElem, K <: RingElem}(a::Integer, b::GenSRPoly{T, K}) = parent(b)(a) - b

-{T <: RingElem, K <: RingElem}(a::fmpz, b::GenSRPoly{T, K}) = parent(b)(a) - b

###############################################################################
#
#   Comparison
#
###############################################################################

function =={T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   if length(a) != length(b)
      return false
   end
   for i = 1:length(a)
      if a.exps[i] != b.exps[i] || a.coeffs[i] != b.coeffs[i]
         return false
      end
   end
   return true
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function =={T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::Int)
   return length(a) == 0 ? b == 0 : a.length == 1 && 
          a.exps[1] == 0 && a.coeffs[1] == b
end

=={T <: RingElem, K <: RingElem}(a::Int, b::GenSRPoly{T, K}) = b == a

function =={T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::fmpz)
   return length(a) == 0 ? b == 0 : a.length == 1 &
          a.exps[1] == 0 && a.coeffs[1] == b
end

=={T <: RingElem, K <: RingElem}(a::fmpz, b::GenSRPoly{T, K}) = b == a

###############################################################################
#
#   Powering
#
###############################################################################

function from_exp(R::Ring, a::UInt)
   return R(fmpz(a))
end

# Implement fps algorithm from "Sparse polynomial powering with heaps" by
# Monagan and Pearce, except that we modify the algorithm to return terms
# in ascending order and we fix some issues in the original algorithm
# http://www.cecm.sfu.ca/CAG/papers/SparsePowering.pdf

function pow_fps{T <: RingElem, K <: RingElem}(f::GenSRPoly{T, K}, k::Int)
   par = parent(f)
   R = base_ring(par)
   m = length(f)
   H = Array(heap_sr, 0) # heap
   I = Array(heap_t, 0) # auxilliary data for heap nodes
   # set up output poly coeffs and exponents (corresponds to h in paper)
   r_alloc = k*(m - 1) + 1
   Rc = Array(T, r_alloc)
   Re = Array(UInt, r_alloc)
   rnext = 1
   # set up g coeffs and exponents (corresponds to g in paper)
   g_alloc = k*(m - 1) + 1
   gc = Array(T, g_alloc)
   ge = Array(UInt, g_alloc)
   gnext = 1
   # set up heap
   gc[1] = f.coeffs[1]^(k-1)
   ge[1] = f.exps[1]*(k - 1)
   Rc[1] = f.coeffs[1]*gc[1]
   Re[1] = f.exps[1]*k
   push!(H, heap_sr(f.exps[2] + ge[1], 1))
   push!(I, heap_t(2, 1, 0))
   Q = Array(Int, 0) # corresponds to Q in paper
   topbit = -1 << (sizeof(Int)*8 - 1)
   mask = ~topbit
   largest = fill(topbit, m) # largest j s.t. (i, j) has been in heap
   largest[2] = 1
   # precompute some values
   fik = Array(T, m)
   for i = 1:m
      fik[i] = from_exp(R, f.exps[i])*(k - 1)
   end
   kp1f1 = k*from_exp(R, f.exps[1])
   gi = Array(T, 1)
   gi[1] = -from_exp(R, ge[1])
   finalexp = f.exps[m]*(k - 1) + f.exps[1]
   # temporary space
   t1 = R()
   C = R() # corresponds to C in paper
   SS = R() # corresponds to S in paper
   temp = R() # temporary space for addmul
   temp2 = R() # temporary space for add
   # begin algorithm
   while !isempty(H)
      exp = H[1].exp
      gnext += 1
      rnext += 1
      if gnext > g_alloc
         g_alloc *= 2
         resize!(gc, g_alloc)
         resize!(ge, g_alloc)
      end
      if rnext > r_alloc
         r_alloc *= 2
         resize!(Rc, r_alloc)
         resize!(Re, r_alloc)
      end
      first = true
      zero!(C) 
      zero!(SS)
      while !isempty(H) && H[1].exp == exp
         x = H[1]
         heappop!(H)
         v = I[x.n]
         largest[v.i] |= topbit
         mul!(t1, f.coeffs[v.i], gc[v.j])
         addeq!(SS, t1)
         if exp <= finalexp
            add!(temp2, fik[v.i], gi[v.j])
            addmul!(C, temp2, t1, temp)
         end
         if first
            ge[gnext] = exp - f.exps[1]
            first = false
         end
         push!(Q, x.n)
         while (xn = v.next) != 0
            v = I[xn]
            largest[v.i] |= topbit
            mul!(t1, f.coeffs[v.i], gc[v.j])
            addeq!(SS, t1)
            if exp <= finalexp
               add!(temp2, fik[v.i], gi[v.j])
               addmul!(C, temp2, t1, temp)
            end
            push!(Q, xn)
         end
      end
      reuse = 0
      while !isempty(Q)
         xn = pop!(Q)
         v = I[xn]
         if v.i < m && largest[v.i + 1] == ((v.j - 1) | topbit)
            I[xn] = heap_t(v.i + 1, v.j, 0)
            heapinsert!(H, I, xn, f.exps[v.i + 1] + ge[v.j]) # either chain or insert v into heap   
            largest[v.i + 1] = v.j
         else
            reuse = xn
         end
         if v.j < gnext - 1 && (largest[v.i] & mask) <  v.j + 1
            if reuse != 0
               I[reuse] = heap_t(v.i, v.j + 1, 0)
               heapinsert!(H, I, reuse, f.exps[v.i] + ge[v.j + 1]) # either chain or insert v into heap
               reuse = 0   
            else
               push!(I, heap_t(v.i, v.j + 1, 0))
               Collections.heappush!(H, heap_sr(f.exps[v.i] + ge[v.j + 1], length(I)))
            end
            largest[v.i] = v.j + 1     
         end
      end
      if C != 0
         temp = divexact(C, from_exp(R, exp) - kp1f1)
         addeq!(SS, temp)
         gc[gnext] = divexact(temp, f.coeffs[1])
         push!(gi, -from_exp(R, ge[gnext]))
         if (largest[2] & topbit) != 0
            push!(I, heap_t(2, gnext, 0))
            Collections.heappush!(H, heap_sr(f.exps[2] + ge[gnext], length(I)))   
            largest[2] = gnext
         end
      end
      if SS != 0
         Rc[rnext] = SS
         Re[rnext] = ge[gnext] + f.exps[1]
         SS = R()
      else
         rnext -= 1
      end
      if C == 0
         gnext -= 1
      end
   end
   resize!(Rc, rnext)
   resize!(Re, rnext)
   return parent(f)(Rc, Re)
end

function ^{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::Int)
   b < 0 && throw(DomainError())
   # special case powers of x for constructing polynomials efficiently
   if length(a) == 0
      return parent(a)()
   elseif length(a) == 1
      return parent(a)([coeff(a, 0)^b], [a.exps[1]*b])
   elseif b == 0
      return parent(a)(1)
   elseif b == 1
      return a
   elseif b == 2
      return a*a
   else
      return pow_fps(a, b)
   end
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divides_monagan_pearce{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   par = parent(a)
   R = base_ring(par)
   m = length(a)
   n = length(b)
   n == 0 && error("Division by zero in divides_monagan_pearce")
   if m == 0
      return true, par()
   end
   H = Array(heap_sr, 0)
   I = Array(heap_t, 0)
   # set up heap
   push!(H, heap_sr(a.exps[1], 1))
   push!(I, heap_t(0, 1, 0))
   q_alloc = max(m - n, n)
   Qc = Array(T, q_alloc)
   Qe = Array(UInt, q_alloc)
   k = 0
   s = n
   c = R()
   qc = R()
   m1 = -R(1)
   mb = -b.coeffs[1]
   Q = Array(Int, 0)
   reuse = Array(Int, 0)
   while !isempty(H)
      exp = H[1].exp
      k += 1
      if k > q_alloc
         q_alloc *= 2
         resize!(Qc, q_alloc)
         resize!(Qe, q_alloc)
      end
      first = true
      @inbounds while !isempty(H) && H[1].exp == exp
         x = H[1]
         heappop!(H)
         v = I[x.n]
         if first
            if exp < b.exps[1]
               return false, par()
            end
            Qe[k] = exp - b.exps[1]
            first = false
         end
         if v.i == 0
            addmul!(qc, a.coeffs[v.j], m1, c)
         else
            addmul!(qc, b.coeffs[v.i], Qc[v.j], c)
         end
         if v.i != 0 || v.j < m
            push!(Q, x.n)
         else
            push!(reuse, x.n)
         end
         while (xn = v.next) != 0
            v = I[xn]
            if v.i == 0
               addmul!(qc, a.coeffs[v.j], m1, c)
            else
               addmul!(qc, b.coeffs[v.i], Qc[v.j], c)
            end
            if v.i != 0 || v.j < m
               push!(Q, xn)
            else
               push!(reuse, xn)
            end
         end
      end
      @inbounds while !isempty(Q)
         xn = pop!(Q)
         v = I[xn]
         if v.i == 0
            I[xn] = heap_t(0, v.j + 1, 0)
            heapinsert!(H, I, xn, a.exps[v.j + 1]) # either chain or insert into heap   
         elseif v.j < k - 1
            I[xn] = heap_t(v.i, v.j + 1, 0)
            heapinsert!(H, I, xn, b.exps[v.i] + Qe[v.j + 1]) # either chain or insert into heap
         elseif v.j == k - 1
            s += 1
            push!(reuse, xn)
         end  
      end
      if qc == 0
         k -= 1
      else
         d2, Qc[k] = divides(qc, mb)
         if !d2
             return false, par()
         end
         for i = 2:s
            if !isempty(reuse)
               xn = pop!(reuse)
               I[xn] = heap_t(i, k, 0)
               heapinsert!(H, I, xn, b.exps[i] + Qe[k]) # either chain or insert into heap
            else
               push!(I, heap_t(i, k, 0))
               Collections.heappush!(H, heap_sr(b.exps[i] + Qe[k], length(I)))
            end                 
         end
         s = 1
      end
      zero!(qc)
   end
   resize!(Qc, k)
   resize!(Qe, k)
   return true, parent(a)(Qc, Qe)
end

function divides{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   d1 = a.exps[a.length]
   d2 = b.exps[b.length] - b.exps[1]
   q_alloc = b.length
   Qe = Array(UInt, q_alloc)
   Qc = Array(T, q_alloc)
   r = a
   k = 0
   while length(r) > 0
      if r.exps[1] + d2 > d1
         return false, parent(a)()
      end
      flag, c = divides(r.coeffs[1], b.coeffs[1])
      if !flag
         return false, parent(a)()
      end
      s = b*c
      if r.exps[1] < b.exps[1]
         return false, parent(a)()
      end
      d = r.exps[1] - b.exps[1]
      k += 1
      if k > q_alloc
         q_alloc *= 2
         resize!(Qe, q_alloc)
         resize!(Qc, q_alloc)
      end
      Qe[k] = d
      Qc[k] = c
      shift_left!(s, d)
      r -= s
   end
   resize!(Qe, k)
   resize!(Qc, k)
   return true, parent(a)(Qc, Qe)
end

function divexact{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   d, q = divides(a, b)
   d == false && error("Not an exact division in divexact")
   return q
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::T)
   len = length(a)
   exps = deepcopy(a.exps)
   coeffs = [divexact(a.coeffs[i], b) for i in 1:len]
   return parent(a)(coeffs, exps)
end

###############################################################################
#
#   Pseudodivision
#
###############################################################################

function pseudodivrem{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   par = parent(a)
   R = base_ring(par)
   m = length(a)
   n = length(b)
   maxn = a.exps[m]
   n == 0 && error("Division by zero in pseudodivrem")
   if a.exps[m] < b.exps[n]
      return par(), a
   end
   H = Array(heap_sr, 0)
   I = Array(heap_t, 0)
   # set up heap
   push!(H, heap_sr(maxn - a.exps[m], 1))
   push!(I, heap_t(0, 1, 0))
   q_alloc = max(m - n + 1, n)
   r_alloc = n
   p_alloc = max(m - n + 1, 1)
   Qc = Array(T, q_alloc)
   Qe = Array(UInt, q_alloc)
   Qpow = Array(Int, q_alloc)
   Rc = Array(T, r_alloc)
   Re = Array(UInt, r_alloc)
   Pc = Array(T, p_alloc)
   Pc[1] = b.coeffs[n]
   p = -1 # current power of lead(b)
   p_max = 1
   k = 0
   l = 0
   s = n
   c = R()
   m1 = -R(1)
   mb = -b.coeffs[n]
   Q = Array(Int, 0)
   reuse = Array(Int, 0)
   while !isempty(H)
      exp = H[1].exp
      qc = R()
      k += 1
      if k > q_alloc
         q_alloc *= 2
         resize!(Qc, q_alloc)
         resize!(Qe, q_alloc)
         resize!(Qpow, q_alloc)
      end
      p += 1
      if p > p_alloc
         p_alloc *= 2
         resize!(Pc, p_alloc)
      end
      if p > p_max
         Pc[p] = Pc[p - 1]*b.coeffs[n]
         p_max = p
      end
      @inbounds while !isempty(H) && H[1].exp == exp
         x = H[1]
         heappop!(H)
         v = I[x.n]
         if v.i == 0
            # cannot use mul! here since Pc[p] is one variable less
            if p > 0
               tc = a.coeffs[m + 1 - v.j]*(-Pc[p])
            else
               tc = -a.coeffs[m + 1 - v.j]
            end
            addeq!(qc, tc)
         else
            if p - Qpow[v.j] - 1 > 0
               tc = Qc[v.j]*Pc[p - Qpow[v.j] - 1]
            else
               tc = Qc[v.j]
            end
            addmul!(qc, b.coeffs[n + 1 - v.i], tc, c)
         end
         if v.i != 0 || v.j < m
            push!(Q, x.n)
         else
            push!(reuse, x.n)
         end
         while (xn = v.next) != 0
            v = I[xn]
            if v.i == 0
               if p > 0
                  tc = a.coeffs[m + 1 - v.j]*(-Pc[p])
               else
                  tc = -a.coeffs[m + 1 - v.j]
               end
               addeq!(qc, tc)
            else
               if p - Qpow[v.j] - 1 > 0
                  tc = Qc[v.j]*Pc[p - Qpow[v.j] - 1]
               else
                  tc = Qc[v.j]
               end
               addmul!(qc, b.coeffs[n + 1 - v.i], tc, c)
            end
            if v.i != 0 || v.j < m
               push!(Q, xn)
            else
               push!(reuse, xn)
            end
         end
      end
      @inbounds while !isempty(Q)
         xn = pop!(Q)
         v = I[xn]
         if v.i == 0
            I[xn] = heap_t(0, v.j + 1, 0)
            heapinsert!(H, I, xn, maxn - a.exps[m - v.j]) # either chain or insert into heap   
         elseif v.j < k - 1
            I[xn] = heap_t(v.i, v.j + 1, 0)
            heapinsert!(H, I, xn, maxn - b.exps[n + 1 - v.i] - Qe[v.j + 1]) # either chain or insert into heap
         elseif v.j == k - 1
            s += 1
            push!(reuse, xn)
         end  
      end
      if qc == 0
         k -= 1
         p -= 1
      else
         d1 = maxn - b.exps[n] >= exp
         texp = maxn - b.exps[n] - exp
         if !d1
            l += 1
            if l >= r_alloc
               r_alloc *= 2
               resize!(Rc, r_alloc)
               resize!(Re, r_alloc)
            end
            Rc[l] = -qc
            Re[l] = maxn - exp
            k -= 1
            p -= 1
         else
            Qc[k] = -qc
            Qe[k] = texp
            Qpow[k] = p
            for i = 2:s
               if !isempty(reuse)
                  xn = pop!(reuse)
                  I[xn] = heap_t(i, k, 0)
                  heapinsert!(H, I, xn, maxn - b.exps[n + 1 - i] - Qe[k]) # either chain or insert into heap
               else
                  push!(I, heap_t(i, k, 0))
                  Collections.heappush!(H, heap_sr(maxn - b.exps[n + 1 - i] - Qe[k], length(I)))
               end
            end                 
            s = 1
         end
      end
   end
   resize!(Qc, k)
   resize!(Qe, k)
   for i = 1:k
      if p > Qpow[i]
         mul!(Qc[i], Qc[i], Pc[p - Qpow[i]])
      else
         Qc[i] = Qc[i]
      end
   end
   resize!(Rc, l)
   resize!(Re, l)
   reverse!(Qc)
   reverse!(Qe)
   reverse!(Rc)
   reverse!(Re)
   if a.exps[m] - b.exps[n] > p
      k1 = reinterpret(Int, a.exps[m] - b.exps[n] - p)
      p1 = b.coeffs[n]^k1
      return parent(a)(Qc, Qe)*p1, parent(a)(Rc, Re)*p1
   else
      return parent(a)(Qc, Qe), parent(a)(Rc, Re)
   end
end

function pseudorem_monagan_pearce{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   par = parent(a)
   R = base_ring(par)
   m = length(a)
   n = length(b)
   maxn = a.exps[m]
   n == 0 && error("Division by zero in pseudodivrem")
   if a.exps[m] < b.exps[n]
      return a
   end
   H = Array(heap_sr, 0)
   I = Array(heap_t, 0)
   # set up heap
   push!(H, heap_sr(maxn - a.exps[m], 1))
   push!(I, heap_t(0, 1, 0))
   q_alloc = max(m - n + 1, n)
   r_alloc = n
   p_alloc = max(m - n + 1, 1)
   Qc = Array(T, q_alloc)
   Qe = Array(UInt, q_alloc)
   Qpow = Array(Int, q_alloc)
   Rc = Array(T, r_alloc)
   Re = Array(UInt, r_alloc)
   Pc = Array(T, p_alloc)
   Pc[1] = b.coeffs[n]
   p = -1 # current power of lead(b)
   p_max = 1
   k = 0
   l = 0
   s = n
   c = R()
   m1 = -R(1)
   mb = -b.coeffs[n]
   Q = Array(Int, 0)
   reuse = Array(Int, 0)
   while !isempty(H)
      exp = H[1].exp
      qc = R()
      k += 1
      if k > q_alloc
         q_alloc *= 2
         resize!(Qc, q_alloc)
         resize!(Qe, q_alloc)
         resize!(Qpow, q_alloc)
      end
      p += 1
      if p > p_alloc
         p_alloc *= 2
         resize!(Pc, p_alloc)
      end
      if p > p_max
         Pc[p] = Pc[p - 1]*b.coeffs[n]
         p_max = p
      end
      @inbounds while !isempty(H) && H[1].exp == exp
         x = H[1]
         heappop!(H)
         v = I[x.n]
         if v.i == 0
            # cannot use mul! here since Pc[p] is one variable less
            if p > 0
               tc = a.coeffs[m + 1 - v.j]*(-Pc[p])
            else
               tc = -a.coeffs[m + 1 - v.j]
            end
            addeq!(qc, tc)
         else
            if p - Qpow[v.j] - 1 > 0
               tc = Qc[v.j]*Pc[p - Qpow[v.j] - 1]
            else
               tc = Qc[v.j]
            end
            addmul!(qc, b.coeffs[n + 1 - v.i], tc, c)
         end
         if v.i != 0 || v.j < m
            push!(Q, x.n)
         else
            push!(reuse, x.n)
         end
         while (xn = v.next) != 0
            v = I[xn]
            if v.i == 0
               if p > 0
                  tc = a.coeffs[m + 1 - v.j]*(-Pc[p])
               else
                  tc = -a.coeffs[m + 1 - v.j]
               end
               addeq!(qc, tc)
            else
               if p - Qpow[v.j] - 1 > 0
                  tc = Qc[v.j]*Pc[p - Qpow[v.j] - 1]
               else
                  tc = Qc[v.j]
               end
               addmul!(qc, b.coeffs[n + 1 - v.i], tc, c)
            end
            if v.i != 0 || v.j < m
               push!(Q, xn)
            else
               push!(reuse, xn)
            end
         end
      end
      @inbounds while !isempty(Q)
         xn = pop!(Q)
         v = I[xn]
         if v.i == 0
            I[xn] = heap_t(0, v.j + 1, 0)
            heapinsert!(H, I, xn, maxn - a.exps[m - v.j]) # either chain or insert into heap   
         elseif v.j < k - 1
            I[xn] = heap_t(v.i, v.j + 1, 0)
            heapinsert!(H, I, xn, maxn - b.exps[n + 1 - v.i] - Qe[v.j + 1]) # either chain or insert into heap
         elseif v.j == k - 1
            s += 1
            push!(reuse, xn)
         end  
      end
      if qc == 0
         k -= 1
         p -= 1
      else
         d1 = maxn - b.exps[n] >= exp
         texp = maxn - b.exps[n] - exp
         if !d1
            l += 1
            if l >= r_alloc
               r_alloc *= 2
               resize!(Rc, r_alloc)
               resize!(Re, r_alloc)
            end
            Rc[l] = -qc
            Re[l] = maxn - exp
            k -= 1
            p -= 1
         else
            Qc[k] = -qc
            Qe[k] = texp
            Qpow[k] = p
            for i = 2:s
               if !isempty(reuse)
                  xn = pop!(reuse)
                  I[xn] = heap_t(i, k, 0)
                  heapinsert!(H, I, xn, maxn - b.exps[n + 1 - i] - Qe[k]) # either chain or insert into heap
               else
                  push!(I, heap_t(i, k, 0))
                  Collections.heappush!(H, heap_sr(maxn - b.exps[n + 1 - i] - Qe[k], length(I)))
               end
            end                 
            s = 1
         end
      end
   end
   resize!(Rc, l)
   resize!(Re, l)
   reverse!(Rc)
   reverse!(Re)
   if a.exps[m] - b.exps[n] > p
      k1 = reinterpret(Int, a.exps[m] - b.exps[n] - p)
      p1 = b.coeffs[n]^k1
      return parent(a)(Rc, Re)*p1
   else
      return parent(a)(Rc, Re)
   end
end

function pseudorem{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   par = parent(a)
   R = base_ring(par)
   m = length(a)
   n = length(b)
   n == 0 && throw(DivideError())
   if a.exps[a.length] < b.exps[b.length]
      return a
   end
   if n > 30
      return pseudorem_monagan_pearce(a, b)
   end
   k = reinterpret(Int, a.exps[a.length] - b.exps[b.length]) + 1
   l = lead(b)
   while a.length > 0 && a.exps[a.length] >= b.exps[b.length]
      s = lead(a)*b
      shift_left!(s, a.exps[a.length] - b.exps[b.length])
      a = a*l - s
      k -= 1
   end
   return a*l^k
end

###############################################################################
#
#   GCD, content and primitive part
#
###############################################################################

function gcd{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   if b.exps[b.length] > a.exps[a.length]
      (a, b) = (b, a)
   end
   if b == 0
      return a
   end
   if b == 1
      return b
   end
   c1 = content(a)
   c2 = content(b)
   c = gcd(c1, c2)
   a = divexact(a, c1)
   b = divexact(b, c2)
   g = one(base_ring(a))
   h = one(base_ring(a))
   while true
      d = reinterpret(Int, a.exps[a.length] - b.exps[b.length])
      r = pseudorem(a, b)
      if r == 0
         break
      end
      if isone(r)
         b = one(parent(a))
         break
      end
      (a, b) = (b, divexact(r, g*h^d))
      g = lead(a)
      if d > 1
         h = divexact(g^d, h^(d - 1))
      else
         h = h^(1 - d)*g^d
      end
   end
   return c*primpart(b)
end

function content{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K})
   z = coeff(a, 0)
   for i = 2:length(a)
      z = gcd(z, coeff(a, i - 1))
   end
   return z
end

function primpart{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K})
   d = content(a)
   return divexact(a, d)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, n::Int)
   if length(a.coeffs) < n
      resize!(a.coeffs, n)
      resize!(a.exps, n)
   end
   return
end

function addmul!{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K}, c::GenSRPoly{T, K}, d::GenSRPoly{T, K})
   t = b*c
   t += a
   a.coeffs = t.coeffs
   a.exps = t.exps
   a.length = t.length
   return
end

function mul!{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K}, c::GenSRPoly{T, K})
   t = b*c
   a.coeffs = t.coeffs
   a.exps = t.exps
   a.length = t.length
   return
end

function addeq!{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K})
   t = a + b
   a.coeffs = t.coeffs
   a.exps = t.exps
   a.length = t.length
   return
end

function add!{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, b::GenSRPoly{T, K}, c::GenSRPoly{T, K})
   t = b + c
   a.coeffs = t.coeffs
   a.exps = t.exps
   a.length = t.length
   return
end

function zero!{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K})
   a.length = 0
   return
end

function shift_left!{T <: RingElem, K <: RingElem}(a::GenSRPoly{T, K}, n::UInt)
   for i = 1:length(a)
      a.exps[i] += n
   end
   return
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

Base.promote_rule{T <: RingElem, K <: RingElem, V <: Integer}(::Type{GenSRPoly{T, K}}, ::Type{V}) = GenSRPoly{T, K}

Base.promote_rule{T <: RingElem, K <: RingElem}(::Type{GenSRPoly{T, K}}, ::Type{T}) = GenSRPoly{T, K}

function promote_rule1{T <: RingElem, K <: RingElem, U <: RingElem}(::Type{GenSRPoly{T, K}}, ::Type{GenSRPoly{U, K}})
   Base.promote_rule(T, GenSRPoly{U, K}) == T ? GenSRPoly{T, K} : Union{}
end

function Base.promote_rule{T <: RingElem, K <: RingElem, U <: RingElem}(::Type{GenSRPoly{T, K}}, ::Type{U})
   Base.promote_rule(T, U) == T ? GenSRPoly{T, K} : promote_rule1(U, GenSRPoly{T, K})
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call{T <: RingElem, K <: RingElem}(a::GenSRPolyRing{T, K}, b::RingElem)
   return a([base_ring(a)(b)], [UInt(0)])
end

function Base.call{T <: RingElem, K <: RingElem}(a::GenSRPolyRing{T, K})
   z = GenSRPoly{T, K}()
   z.parent = a
   return z
end

function Base.call{T <: RingElem, K <: RingElem}(a::GenSRPolyRing{T, K}, b::Integer)
   z = GenSRPoly{T, K}(base_ring(a)(b))
   z.parent = a
   return z
end

function Base.call{T <: RingElem, K <: RingElem}(a::GenSRPolyRing{T, K}, b::T)
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = GenSRPoly{T, K}(b)
   z.parent = a
   return z
end

function Base.call{T <: RingElem, K <: RingElem}(a::GenSRPolyRing{T, K}, b::GenSRPoly{T, K})
   parent(b) != a && error("Unable to coerce polynomial")
   return b
end

function Base.call{T <: RingElem, K <: RingElem}(a::GenSRPolyRing{T, K}, b::Array{T, 1}, m::Array{UInt, 1})
   if length(b) > 0
      parent(b[1]) != base_ring(a) && error("Unable to coerce to polynomial")
   end
   z = GenSRPoly{T, K}(b, m)
   z.parent = a
   return z
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

doc"""
    SRPolynomialRing(R::Ring, s::Array{String, 1}; cached::Bool = true)
> Given a base ring `R` and an array of strings `s` specifying how the
> generators (variables) should be printed, return a tuple `S, x1, x2, ...`
> representing the new polynomial ring $T = R[x1, x2, ...]$ and the generators
> $x1, x2, ...$ of the polynomial ring. By default the parent object `T` will
> depend only on `R` and `x1, x2, ...` and will be cached. Setting the optional
> argument `cached` to `false` will prevent the parent object `T` from being
> cached.
"""
function SRPolynomialRing(R::Ring, s::Array{ASCIIString{}, 1}; cached::Bool = true)
   U = [Symbol(s[i]) for i in length(s):-1:1]
   T = elem_type(R)

   R1 = R
   T1 = elem_type(R)
   K = T1
   for i = 1:length(U)
      R1 = GenSRPolyRing{T1, K}(R1, U[i], i, cached)
      T1 = elem_type(R1)
   end

   return tuple(R1, gens(R1)...)
end
