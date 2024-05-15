@testset "WeakKeyIdDict" begin
  import AbstractAlgebra: WeakKeyIdDict
  A = [1]
  B = [2]
  C = [3]

  # construction
  wkd = WeakKeyIdDict()
  wkd[A] = 2
  wkd[B] = 3
  wkd[C] = 4
  dd = convert(Dict{Any,Any},wkd)
  @test WeakKeyIdDict(dd) == wkd
  @test convert(WeakKeyIdDict{Any, Any}, dd) == wkd
  @test isa(WeakKeyIdDict(dd), WeakKeyIdDict{Any,Any})

  # test many constructors without type parameters specified
  @test WeakKeyIdDict(A=>2, B=>3, C=>4) == wkd
  @test isa(WeakKeyIdDict(A=>2, B=>3, C=>4), WeakKeyIdDict{Vector{Int},Int})
  @test WeakKeyIdDict(a=>i+1 for (i,a) in enumerate([A,B,C]) ) == wkd
  @test WeakKeyIdDict([(A,2), (B,3), (C,4)]) == wkd
  @test WeakKeyIdDict(Pair(A,2), Pair(B,3), Pair(C,4)) == wkd

  # test many constructors with type parameters specified
  @test WeakKeyIdDict{Vector{Int},Int}(A=>2, B=>3, C=>4) == wkd
  @test isa(WeakKeyIdDict{Vector{Int},Int}(A=>2, B=>3, C=>4), WeakKeyIdDict{Vector{Int},Int})
  @test WeakKeyIdDict{Vector{Int},Int}(a=>i+1 for (i,a) in enumerate([A,B,C]) ) == wkd
  @test WeakKeyIdDict{Vector{Int},Int}([(A,2), (B,3), (C,4)]) == wkd
  @test WeakKeyIdDict{Vector{Int},Int}(Pair(A,2), Pair(B,3), Pair(C,4)) == wkd

  # test more constructors with mixed types
  @test isa(WeakKeyIdDict(A=>2, B=>3, C=>"4"), WeakKeyIdDict{Vector{Int},Any})
  @test isa(WeakKeyIdDict(A=>2, B=>3, "C"=>4), WeakKeyIdDict{Any,Int})
  @test isa(WeakKeyIdDict(A=>2, B=>3, "C"=>"4"), WeakKeyIdDict{Any,Any})

  @test copy(wkd) == wkd

  @test length(wkd) == 3
  @test !isempty(wkd)
  @test haskey(wkd, C)
  @test getkey(wkd, C, 123) === C
  res = pop!(wkd, C)
  @test res == 4
  @test C ∉ keys(wkd)
  @test 4 ∉ values(wkd)
  @test !haskey(wkd, C)
  @test length(wkd) == 2
  @test !isempty(wkd)
  @test 47 == pop!(wkd, C, 47)
  @test getkey(wkd, C, 123) == 123
  wkd = filter!( p -> p.first != B, wkd)
  @test B ∉ keys(wkd)
  @test 3 ∉ values(wkd)
  @test length(wkd) == 1
  @test WeakKeyIdDict(Pair(A, 2)) == wkd
  @test !isempty(wkd)

  wkd[A] = 42
  @test wkd[A] == 42

  wkd = WeakKeyIdDict(A=>2, B=>3, C=>4)
  map!(x -> x + 1, values(wkd))
  @test WeakKeyIdDict(A=>3, B=>4, C=>5) == wkd

  wkd = WeakKeyIdDict(A=>2, B=>3, C=>4)
  @test delete!(wkd, A) == WeakKeyIdDict(B=>3, C=>4)
  @test delete!(wkd, A) == WeakKeyIdDict(B=>3, C=>4)  # deleting the same key twice works
  @test delete!(wkd, C) == WeakKeyIdDict(B=>3)
  @test delete!(wkd, B) == WeakKeyIdDict()
  # adding stuff back is OK
  wkd[A] = 2
  wkd[B] = 3
  wkd[C] = 4
  @test wkd == WeakKeyIdDict(A=>2, B=>3, C=>4)

  wkd = WeakKeyIdDict(A=>2)
  @test get(wkd, A, 17) == 2
  @test get!(wkd, A, 17) == 2
  @test get(wkd, B, 17) == 17
  @test length(wkd) == 1
  @test get!(wkd, B, 17) == 17
  @test length(wkd) == 2

  wkd = WeakKeyIdDict(A=>2)
  @test get(() -> 23, wkd, A) == 2
  @test get!(() -> 23, wkd, A) == 2
  @test get(() -> 23, wkd, B) == 23
  @test length(wkd) == 1
  @test get!(() -> 23, wkd, B) == 23
  @test length(wkd) == 2

  wkd = empty!(wkd)
  @test wkd == empty(wkd)
  @test wkd == empty(wkd)
  @test typeof(wkd) == typeof(empty(wkd))
  @test length(wkd) == 0
  @test isempty(wkd)
  @test isa(wkd, WeakKeyIdDict)
  @test WeakKeyIdDict() == WeakKeyIdDict(())

  # bad iterable argument
  @test_throws ArgumentError WeakKeyIdDict([1, 2, 3])

  # integers can't be arguments
  @test_throws ErrorException WeakKeyIdDict([1=>2])

  # WeakKeyIdDict does not convert keys
  @test_throws ArgumentError WeakKeyIdDict{Int,Any}(5.0=>1)

  # iterator
  wkd = WeakKeyIdDict(A=>2, B=>3, C=>4)
  @test Set(collect(wkd)) == Set([A=>2, B=>3, C=>4])
  @test 2+3+4 == sum(v for (k,v) in wkd)

  # WeakKeyIdDict hashes with object-id
  AA = copy(A)
  GC.@preserve A AA begin
    wkd = WeakKeyIdDict(A=>1, AA=>2)
    @test length(wkd)==2
    kk = collect(keys(wkd))
    @test kk[1]==kk[2]
    @test kk[1]!==kk[2]
  end

  # WeakKeyIdDict compares to other dicts:
  @test IdDict(A=>1)!=WeakKeyIdDict(A=>1)
  @test Dict(A=>1)==WeakKeyIdDict(A=>1)
  @test Dict(copy(A)=>1)!=WeakKeyIdDict(A=>1)

end

@testset "WeakKeyIdDict.lock" begin
  A = [1]
  B = [2]
  C = [3]
  wkd = WeakKeyIdDict(A=>2, B=>3, C=>4)
  @test !islocked(wkd)
  lock(wkd)
  @test islocked(wkd)
  unlock(wkd)
  @test !islocked(wkd)
end

@testset "WeakKeyIdDict.gc" begin

  # issue #26939
  d26939 = WeakKeyIdDict()
  d26939[big"1.0" + 1.1] = 1
  GC.gc() # make sure this doesn't segfault

  wkd = WeakKeyIdDict([42]=>2, [43]=>3, [44]=>4)
  for k in keys(wkd)
    delete!(wkd, k)
  end
  @test isempty(wkd)
  GC.gc()  # try to get it to evict some weak references
  @test isempty(wkd)
end


# verify that garbage collection takes care of our weak references
#
# This test somehow doesn't work if it is inside a @testset as then the
# WeakRef is not collected for some reason I don't quite understand at
# this point... :-(
_tmp_key = [1]
wkd = WeakKeyIdDict(_tmp_key => 1)
let tmp = [ 42 ]
  @test length(wkd) == 1
  wkd[tmp] = 2
  @test length(wkd) == 2
end
# at this point there is no strong reference left to the vector [42]
# previously reachable via tmp
GC.gc(true)

@test length(wkd) == 1
@test length(keys(wkd)) == 1
@test WeakKeyIdDict(_tmp_key => 1) == wkd
