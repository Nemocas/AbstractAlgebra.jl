###############################################################################
#
#  Conformance tests for the group and group element interfaces.
#
#
#  The code in this file is based on v0.4 of the GroupsCore.jl package, which
#  is copyright (c) 2021 Marek Kaluba <kalmar@amu.edu.pl> and contributors,
#  and licensed under the MIT license <https://spdx.org/licenses/MIT.html>.
#
#  The interfaces have since then diverged, with several breaking changes.
#  Function and type names differ (e.g. `istrivial` vs. `is_trivial`), and
#  also semantics (e.g. definition of equality for group elements in
#  GroupsCore allows returning false for mathematically equal elements when it
#  is "too hard" to prove equality; we don't).
#
###############################################################################


function test_Group_interface(G::Group)
  @testset "Group interface" begin
#        @testset "Iteration protocol" begin
#            IS = Base.IteratorSize(typeof(G))
#            if IS isa Base.IsInfinite
#                @test is_finite(G) == false
#            else
#                isfiniteG = false
#                if IS isa Base.HasLength || IS isa Base.HasShape
#                    @test is_finite(G) == true
#                    isfiniteG = true
#                else
#                    @test IS isa Base.SizeUnknown
#                    try
#                        @test is_finite(G) isa Bool
#                        isfiniteG = is_finite(G)
#                    catch e
#                        @test e isa AbstractAlgebra.InfiniteOrderError
#                        isfiniteG = false
#                    end
#                end
#
#                if isfiniteG
#                    @test length(G) isa Int
#                    @test length(G) > 0
#
#                    @test elem_type(G) <: GroupElem
#                    @test one(G) isa elem_type(G)
#
#                    if has_gens(G)
#                        @test first(iterate(G)) isa elem_type(G)
#                        _, s = iterate(G)
#                        @test first(iterate(G, s)) isa elem_type(G)
#                        @test isone(first(G))
#                    end
#                end
#            end
#        end

      @testset "Group generators" begin
          @test has_gens(G) isa Bool

          if has_gens(G)
              @test ngens(G) isa Int
              @test gens(G) isa AbstractVector{elem_type(G)}
              @test length(gens(G)) == ngens(G)
              if ngens(G) > 0
                  @test first(gens(G)) == gen(G, 1)
                  @test last(gens(G)) == gen(G, ngens(G))
              end
          end
      end

      @testset "order, rand" begin
          if is_finite(G)
              ord = order(BigInt, G)
              @test ord isa BigInt
              if ord < typemax(Int16)
                  @test order(Int16, G) isa Int16
              else
                  @test_throws InexactError order(Int16, G)
              end
              @test order(G) >= 1
              @test is_trivial(G) == (order(G) == 1)
          else
              @test_throws AbstractAlgebra.InfiniteOrderError order(G)
              @test !is_trivial(G)
          end

          @test rand(G) isa GroupElem
          @test rand(G, 2) isa AbstractVector{elem_type(G)}
          g, h = rand(G, 2)
          @test parent(g) === parent(h) === G

#            @test GroupsCore.rand_pseudo(G) isa elem_type(G)
#            @test GroupsCore.rand_pseudo(G, 2, 2) isa AbstractMatrix{elem_type(G)}
#
#            g, h = GroupsCore.rand_pseudo(G, 2)
#            @test parent(g) === parent(h) === G
      end
  end
end

function test_GroupElem_interface(g::GEl, h::GEl) where {GEl<:GroupElem}

  @testset "GroupElem interface" begin

      @testset "Parent methods" begin
          @test parent(g) isa Group
          @test parent(g) === parent(h)
          G = parent(g)

          @test elem_type(G) == typeof(g)

          @test one(g) isa elem_type(G)

          @test one(G) == one(g) == one(h)

          @test isone(one(G))
      end

      @testset "Equality, deepcopy && hash" begin
          @test (g == h) isa Bool
          @test isequal(g, h) isa Bool

          @test g == g
          @test isequal(g, g)

          if g != h
              @test !isequal(g, h)
          end

          @test deepcopy(g) isa typeof(g)
          @test deepcopy(g) == g
          k = deepcopy(g)
          @test parent(k) === parent(g)
          @test hash(g) isa UInt
          @test hash(g) == hash(k)

          if isequal(g, h)
              @test hash(g) == hash(h)
          end
      end

      @testset "Group operations" begin
          old_g, old_h = deepcopy(g), deepcopy(h)

          # check that the default operations don't mutate their arguments
          @test inv(g) isa typeof(g)
          @test (g, h) == (old_g, old_h)

          @test g * h isa typeof(g)
          @test (g, h) == (old_g, old_h)

          @test g^2 == g * g
          @test (g, h) == (old_g, old_h)

          @test g^-3 == inv(g) * inv(g) * inv(g)
          @test (g, h) == (old_g, old_h)

          @test (g * h)^-1 == inv(h) * inv(g)
          @test (g, h) == (old_g, old_h)

          @test conj(g, h) == inv(h) * g * h
          @test (g, h) == (old_g, old_h)

          @test ^(g, h) == inv(h) * g * h
          @test (g, h) == (old_g, old_h)

          @test comm(g, h) == g^-1 * h^-1 * g * h
          @test (g, h) == (old_g, old_h)

          @test comm(g, h, g) == conj(inv(g), h) * conj(conj(g, h), g)
          @test (g, h) == (old_g, old_h)

          @test isone(g * inv(g)) && isone(inv(g) * g)
          @test (g, h) == (old_g, old_h)

          @test g / h == g * inv(h)
          @test (g, h) == (old_g, old_h)

          @test div_right(g, h) == g * inv(h)
          @test (g, h) == (old_g, old_h)

          @test g \ h == inv(g) * h
          @test (g, h) == (old_g, old_h)

          @test div_left(g, h) == inv(h) * g
          @test (g, h) == (old_g, old_h)
      end

      @testset "Misc GroupElem methods" begin
          @test one(g) isa typeof(g)
          @test isone(g) isa Bool
          @test isone(one(g))

          @test is_finite_order(g) isa Bool

          if is_finite_order(g)
              @test order(Int16, g) isa Int16
              @test order(BigInt, g) isa BigInt
              @test order(g) >= 1
              if is_finite(parent(g))
                  @test iszero(order(parent(g)) % order(g))
              end
              if !isone(g) && !isone(g^2)
                  @test order(g) > 2
              end
              @test order(inv(g)) == order(g)
              @test order(one(g)) == 1
          else
              @test_throws AbstractAlgebra.InfiniteOrderError order(g)
          end

          @test similar(g) isa typeof(g)
      end

      @testset "In-place operations" begin
          old_g, old_h = deepcopy(g), deepcopy(h)
          out = similar(g)

          @test isone(one!(g))
          g = deepcopy(old_g)

          @test inv!(out, g) == inv(old_g)
          @test g == old_g
          @test inv!(out, g) == inv(old_g)
          g = deepcopy(old_g)

          @testset "mul!" begin
              @test mul!(out, g, h) == old_g * old_h
              @test (g, h) == (old_g, old_h)

              @test mul!(out, g, h) == old_g * old_h
              @test (g, h) == (old_g, old_h)

              @test mul!(g, g, h) == old_g * old_h
              @test h == old_h
              g = deepcopy(old_g)

              @test mul!(h, g, h) == old_g * old_h
              @test g == old_g
              h = deepcopy(old_h)

              @test mul!(g, g, g) == old_g * old_g
              g = deepcopy(old_g)
          end

          @testset "conj!" begin
              res = old_h^-1 * old_g * old_h
              @test conj!(out, g, h) == res
              @test (g, h) == (old_g, old_h)

              @test conj!(g, g, h) == res
              @test h == old_h
              g = deepcopy(old_g)

              @test conj!(h, g, h) == res
              @test g == old_g
              h = deepcopy(old_h)

              @test conj!(g, g, g) == old_g
              g = deepcopy(old_g)
          end

          @testset "comm!" begin
              res = old_g^-1 * old_h^-1 * old_g * old_h

              @test comm!(out, g, h) == res
              @test (g, h) == (old_g, old_h)

              @test comm!(out, g, h) == res
              @test (g, h) == (old_g, old_h)

              @test comm!(g, g, h) == res
              @test h == old_h
              g = deepcopy(old_g)

              @test comm!(h, g, h) == res
              @test g == old_g
              h = deepcopy(old_h)
          end

          @testset "div_[left|right]!" begin
              res = g * h^-1
              @test div_right!(out, g, h) == res
              @test (g, h) == (old_g, old_h)

              @test div_right!(g, g, h) == res
              @test h == old_h
              g = deepcopy(old_g)

              @test div_right!(h, g, h) == res
              @test g == old_g
              h = deepcopy(old_h)

              @test div_right!(g, g, g) == one(g)
              g = deepcopy(old_g)


              res = h^-1 * g
              @test div_left!(out, g, h) == res
              @test (g, h) == (old_g, old_h)

              @test div_left!(g, g, h) == res
              @test h == old_h
              g = deepcopy(old_g)

              @test div_left!(h, g, h) == res
              @test g == old_g
              h = deepcopy(old_h)

              @test div_left!(g, g, g) == one(g)
              g = deepcopy(old_g)
          end
      end

      @testset "Broadcasting" begin
          @test g .* [g, h] == [g * g, g * h]
          G = parent(g)
          if has_gens(G)
              @test g .* gens(G) == [g * x for x in gens(G)]
          end
      end
  end
end
