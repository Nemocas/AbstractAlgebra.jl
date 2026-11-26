using Test

@testset "ProbablyCorrect foundations" begin
  add_dodgy_func(:JustForTesting)
  set_dodgy_mode(true)
  set_dodgy_mode(:JustForTesting, false)
  @test !get_dodgy_mode(:JustForTesting)
  set_dodgy_mode(false)
  @test !get_dodgy_mode(:JustForTesting)
  function test_func(n::Int)
    (n >= 0) && return n
    get_dodgy_mode(:JustForTesting) && @record_dodgy_step(:test_func, (n,), -n)
    return -n
  end
  @test is_empty(dodgy_steps_get())
  test_func(1);
  @test is_empty(dodgy_steps_get())
  test_func(-1);
  @test is_empty(dodgy_steps_get())
  
  set_dodgy_mode(:JustForTesting)
  @test get_active_dodgy_funcs() == [:JustForTesting]
  test_func(1);
  @test is_empty(dodgy_steps_get())
  test_func(-1);
  @test length(dodgy_steps_get()) == 1

  dodgy_steps_clear()
  @test is_empty(dodgy_steps_get())

  set_dodgy_mode(:JustForTesting, false)
  @test is_empty(dodgy_steps_get())
  test_func(1);
  @test is_empty(dodgy_steps_get())
  test_func(-1);
  @test is_empty(dodgy_steps_get())
  
  set_dodgy_mode()
  test_func(1);
  @test is_empty(dodgy_steps_get())
  test_func(-1);
  @test length(dodgy_steps_get()) == 1

  dodgy_steps_clear()
  set_dodgy_mode(false)
end
