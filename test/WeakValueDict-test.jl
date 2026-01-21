WeakValueDict = AbstractAlgebra.WeakValueDict

function test_weak_cache(T, kreps, ireps)
   d = T{BigInt, BigInt}()
   # keys missing from/present in ddeg are definitely missing from/present in d
   ddef  = Dict{BigInt, BigInt}()
   # if a key is in d, its value should match that of dmay
   dmay  = Dict{BigInt, BigInt}()
   for k in 1:kreps
      for i in 1:ireps
         for j in 1:2
            x = BigInt(rand(1:999))
            y = BigInt(rand(1:999))
            ddef[x] = y
            d[x] = y
            y = deepcopy(y)
            @test d[x] == y
            @test ddef[x] == y
         end
         x = BigInt(rand(1:999))
         delete!(ddef, x)
         delete!(d, x)
         @test !haskey(d, x)
         x = BigInt(rand(1:999))
         @test get(d, x, BigInt(0)) == get(ddef, x, BigInt(0))

         @test length(d) >= length(ddef)

         for j in 1:2
            x = BigInt(-rand(1:999))
            y = BigInt(-rand(1:999))
            dmay[x] = y
            d[x] = deepcopy(y)
            @test get(d, x, y) == dmay[x]
         end
         x = BigInt(-rand(1:999))
         delete!(dmay, x)
         delete!(d, x)
         @test !haskey(d, x)
         x = BigInt(rand(1:999))
         y = get(dmay, x, nothing)
         @test y == nothing || get(d, x, y) == y
      end
      GC.gc(true)
   end

   empty!(d)
   empty!(ddef)
   empty!(dmay)
   for k in 1:kreps
      for i in 1:ireps
         for j in 1:2
            x = BigInt(rand(1:999))
            y = BigInt(rand(1:999))
            ddef[x] = y
            delete!(d, x)
            @test (get!(d, x) do; return y; end) == y
            y = deepcopy(y)
            @test d[x] == y
            @test ddef[x] == y
         end
         x = BigInt(rand(1:999))
         delete!(ddef, x)
         delete!(d, x)
         @test !haskey(d, x)
         x = BigInt(rand(1:999))
         @test get(d, x, BigInt(0)) == get(ddef, x, BigInt(0))

         @test length(d) >= length(ddef)

         for j in 1:2
            x = BigInt(-rand(1:999))
            y = BigInt(-rand(1:999))
            dmay[x] = y
            delete!(d, x)
            @test (get!(d, x) do; return deepcopy(y); end) == dmay[x]
            @test (get!(d, Int(x)) do; return deepcopy(y); end) == dmay[x]
         end
         x = BigInt(-rand(1:999))
         delete!(dmay, x)
         delete!(d, x)
         @test !haskey(d, x)
         x = BigInt(rand(1:999))
         y = get(dmay, x, nothing)
         @test y == nothing || get(d, x, y) == y
      end
      GC.gc(true)
   end

   empty!(d)
   if T === AbstractAlgebra.WeakValueCache
      AbstractAlgebra.rehash!(d, 5) # test internal function
   end
   empty!(ddef)
   x = BigInt(1)
   y = BigInt(2)
   GC.@preserve x y begin
      ddef[1] = x
      ddef[2] = y

      # test get! where default modifies d
      get!(d, 2) do; d[1] = x; return y; end
      @test d[1] == x
      @test d[2] == y
      @test length(string(d)) > 3

      @test_throws KeyError d[3]
   end

   d[3] = BigInt(4)
   GC.gc(true)
   try z = d[3]
      @test z == 4
   catch e
      @test e isa KeyError
   end

   empty!(d)
   x = BigInt(1)
   y = BigInt(2)
   GC.@preserve x y begin
      d[x] = y
      @test (get(d, x) do; return BigInt(3); end) isa BigInt

      @test pop!(d, x) == y
      @test_throws KeyError pop!(d, x)
      @test pop!(d, x, y) === y
   end

   x = BigInt(1)
   y = BigInt(2)
   z = BigInt(3)
   GC.@preserve x y z begin
      D = T{BigInt, BigInt}(1=>x)
      @test D[1] === x
      @test_throws KeyError D[2]

      D = T{BigInt, BigInt}(1=>x, 2=>y, 3=>z)
      @test D[1] === x
      @test D[2] === y
      @test D[3] === z
      @test_throws KeyError D[4]

      D = T{BigInt, BigInt}([(1,x), (2,y), (3,z)])
      @test D[1] === x
      @test D[2] === y
      @test D[3] === z
      @test_throws KeyError D[4]
   end
end


@testset "WeakValueCache" begin
   test_weak_cache(AbstractAlgebra.WeakValueCache, 30, 20)
end

@testset "WeakValueDict" begin
    test_weak_cache(AbstractAlgebra.WeakValueDict, 30, 20)

    A = [1]
    B = [2]
    C = [3]

    # construction
    wkd = WeakValueDict()
    wkd[2] = A
    wkd[3] = B
    wkd[4] = C
    dd = convert(Dict{Any,Any},wkd)
    @test WeakValueDict(dd) == wkd
    @test convert(WeakValueDict{Any, Any}, dd) == wkd
    @test isa(WeakValueDict(dd), WeakValueDict{Any,Any})
    @test WeakValueDict(2=>A, 3=>B, 4=>C) == wkd
    @test isa(WeakValueDict(2=>A, 3=>B, 4=>C), WeakValueDict{Int, Vector{Int}})
    @test WeakValueDict((i+1)=>a for (i,a) in enumerate([A,B,C]) ) == wkd
    @test WeakValueDict([(2,A), (3,B), (4,C)]) == wkd
    @test WeakValueDict{Int, typeof(A)}(Pair(2,A), Pair(3,B), Pair(4,C)) == wkd
    @test WeakValueDict(Pair(2,A), Pair(3,B), Pair(4,C)) == wkd
    D = [[4.0]]
    @test WeakValueDict(Pair(2,A), Pair(3,B), Pair(4.0,D)) isa WeakValueDict{Any, Any}
    @test isa(WeakValueDict(Pair(2,A), Pair(3.0,B), Pair(4,C)), WeakValueDict{Any, Vector{Int}})
    @test copy(wkd) == wkd

    @test length(wkd) == 3
    @test !isempty(wkd)
    res = pop!(wkd, 3)
    @test res == B
    @test length(wkd) == 2
    res = pop!(wkd, 3, C)
    @test res == C
    @test 3 ∉ keys(wkd)
    @test B ∉ values(wkd)
    @test length(wkd) == 2
    @test !isempty(wkd)
    wkd = filter!( p -> p.first != 2, wkd)
    @test 3 ∉ keys(wkd)
    @test B ∉ values(wkd)
    @test length(wkd) == 1
    @test WeakValueDict(Pair(4, C)) == wkd
    @test !isempty(wkd)

    wkd = empty!(wkd)
    @test wkd == empty(wkd)
    @test typeof(wkd) == typeof(empty(wkd))
    @test length(wkd) == 0
    @test isempty(wkd)
    @test isa(wkd, WeakValueDict)

    @test_throws ArgumentError WeakValueDict([1, 2, 3])

    wkd = WeakValueDict(1=>A)
    @test delete!(wkd, 1) == empty(wkd)
    @test delete!(wkd, 1) === wkd

    GC.@preserve A B C D nothing
end

