using Pkg

if VERSION>=v"1.9"
  @testset "Banners" begin

    function run_repl_code(code::String, proj::String)
      bin = Base.julia_cmd()
      opts = ["--project=$proj", "-i", "-e", "$code; exit();"]
      cmd = Cmd(`$bin $opts`, ignorestatus=true)
      outs = IOBuffer()
      errs = IOBuffer()
      proc = run(pipeline(`$cmd`, stderr=errs, stdout=outs))
      result = String(take!(outs))
      err = String(take!(errs))
      return result, err, proc.exitcode
    end

    # Set up a separate temporary project for some modules that depend on each
    # other with some of them showing banners. Chain is
    # AA -> ModA -> ModB -> ModC
    path = dirname(@__FILE__)
    aadir = dirname(dirname(pathof(AbstractAlgebra)))
    modadir = joinpath(path, "ModA")
    modbdir = joinpath(path, "ModB")
    modcdir = joinpath(path, "ModC")

    # generate temp project
    td = mktempdir()
    code = """
    using Pkg;
    Pkg.develop(path=raw"$aadir");
    Pkg.develop(path=raw"$modadir");
    Pkg.develop(path=raw"$modbdir");
    Pkg.develop(path=raw"$modcdir");
    Pkg.precompile();
    """
    out,err,exitcode = run_repl_code(code, td)
    res = @test exitcode == 0
    if res isa Test.Fail
      println("out\n$out")
      println("err\n$err")
    end

    # Banner of ModA shows
    out, err = run_repl_code("using ModA;", td)
    res = @test strip(out) == "Banner of ModA"
    if res isa Test.Fail
      println("out\n$out")
      println("err\n$err")
    end

    # Banner of ModB shows, but ModA is supressed
    out, err = run_repl_code("using ModB;", td)
    res = @test strip(out) == "Banner of ModB" broken=VERSION>=v"1.11.0-DEV"
    if res isa Test.Fail
      println("out\n$out")
      println("err\n$err")
    end

    # Banner of ModB shows, but ModA is supressed, even if ModA is specifically
    # used after ModB
    out, err = run_repl_code("using ModB; using ModA;", td)
    res = @test strip(out) == "Banner of ModB" broken=VERSION>=v"1.11.0-DEV"
    if res isa Test.Fail
      println("out\n$out")
      println("err\n$err")
    end

    # Banner does not show when our module is a dependency
    out, err = run_repl_code("using ModC;", td)
    res = @test strip(out) == "" broken=VERSION>=v"1.11.0-DEV"
    if res isa Test.Fail
      println("out\n$out")
      println("err\n$err")
    end
  end
end
