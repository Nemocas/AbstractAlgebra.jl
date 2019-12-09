
@testset "Generic.Mat rectangular solve tests." begin

    @testset "Generic.Mat.solve_lu..." begin
        S = QQ

        @testset "Consistent square solve tests." begin            
            for dim = 1:5
                R = MatrixSpace(S, dim, dim)
                U = MatrixSpace(S, dim, rand(1:5))

                M = randmat_with_rank(R, dim, -100:100)
                b = rand(U, -100:100)

                x = Generic.solve_lu(M, b)

                @test M*x == b
            end
        end

        @testset "Consistent rectangular solve tests." begin
            for rdim = 1:3
                for cdim = 1:3
                    R = MatrixSpace(S, rdim, cdim)
                    U = MatrixSpace(S, cdim, rand(1:5))
                    
                    M = randmat_with_rank(R, min(rdim,cdim), -10:10)
                    x = rand(U, -10:10)

                    b = M*x

                    y = Generic.solve_lu(M, b)

                    @test M*y == b
                end
            end

        end
        
        @testset "Inconsistent solve error tests." begin
            for rdim = 1:3
                for cdim = 1:3
                    Mn = MatrixSpace(S, rdim, rdim)
                    Mm = MatrixSpace(S, cdim, cdim)

                    R = MatrixSpace(S, rdim, cdim)

                    t = rand(1:5)
                    U = MatrixSpace(S, rdim, t)
                    
                    # Assgin a random test matrix.
                    D = R()
                    for j = 1:min(rdim, cdim)-1
                        D[j,j] = rand(-100:100)
                    end

                    # Choose vector space automorphisms.
                    g = randmat_with_rank(Mn, rdim, -10:10)
                    h = randmat_with_rank(Mm, cdim, -10:10)

                    # Note no solution to `Dx=b` is possible, as `b` at least one column of
                    # `b` has a non-constant final coordinate.
                    bbad = rand(U, -10:10); bbad[rdim,rand(1:t)] = 1

                    #Act!
                    M = g*D*h
                    bbad = g*bbad

                    #@test Generic.solve_lu(M,bbad)
                    @test_throws DomainError Generic.solve_lu(M,bbad)
                end
            end
        end
        
        @testset "Solves with zero matrix inputs." begin

            A = matrix(QQ, fill(QQ(0),4,4))
            b = matrix(QQ, hcat([[rand(-10:10) for j=1:4] for i=1:1]...))
            b[1,1] = 1
            
            @test_throws DomainError Generic.solve_lu(A,b)

            R = MatrixSpace(S, 5, 4)
            A = rand(R, -10:10)
            b = matrix(QQ, fill(0,5,4))
            x = Generic.solve_lu(A,b)
            @test A*x == b

            A = matrix(QQ, fill(QQ(0),5,4))
            b = matrix(QQ, fill(0,5,4))
            x = Generic.solve_lu(A,b)
            @test A*x == b
        end
        
        
        S, y = PolynomialRing(ZZ, "y")
        K = FractionField(S)

        for dim = 1:5
            R = MatrixSpace(S, dim, dim)
            U = MatrixSpace(S, dim, rand(1:5))

            M = randmat_with_rank(R, dim, 0:5, -100:100)
            b = rand(U, 0:5, -100:100);

            MK = matrix(K, elem_type(K)[ K(M[i, j]) for i in 1:nrows(M), j in 1:ncols(M) ])
            bK = matrix(K, elem_type(K)[ K(b[i, j]) for i in 1:nrows(b), j in 1:ncols(b) ])

            x = Generic.solve_lu(MK, bK)

            @test MK*x == bK
        end
        
    end

    @testset "Generic.Mat.solve_fflu..." begin
        S = ZZ

        @testset "Consistent square solve tests." begin            
            for dim = 1:7
                R = MatrixSpace(S, dim, dim)
                U = MatrixSpace(S, dim, rand(1:5))

                M = randmat_with_rank(R, dim, -100:100)
                b = rand(U, -100:100)

                x,d = Generic.solve_fflu(M, b)

                @test M*x == d*b
            end
        end

        
        @testset "Consistent rectangular solve tests." begin
            for rdim = 1:3
                for cdim = 1:3
                    R = MatrixSpace(S, rdim, cdim)
                    U = MatrixSpace(S, cdim, rand(1:5))
                    
                    M = randmat_with_rank(R, min(rdim,cdim), -10:10)
                    x = rand(U, -10:10)

                    b = M*x

                    y,d = Generic.solve_fflu(M, b)

                    @test M*y == d*b
                end
            end

        end

        
        @testset "Inconsistent solve error tests." begin
            for rdim = 1:3
                for cdim = 1:3
                    Mn = MatrixSpace(S, rdim, rdim)
                    Mm = MatrixSpace(S, cdim, cdim)

                    R = MatrixSpace(S, rdim, cdim)

                    t = rand(1:5)
                    U = MatrixSpace(S, rdim, t)
                    
                    # Assgin a random test matrix.
                    D = R()
                    for j = 1:min(rdim, cdim)-1
                        D[j,j] = rand(-100:100)
                    end

                    # Choose vector space automorphisms.
                    g = randmat_with_rank(Mn, rdim, -10:10)
                    h = randmat_with_rank(Mm, cdim, -10:10)

                    # Note no solution to `Dx=b` is possible, as `b` at least one column of
                    # `b` has a non-constant final coordinate.
                    bbad = rand(U, -10:10); bbad[rdim,rand(1:t)] = 1

                    #Act!
                    M = g*D*h
                    bbad = g*bbad

                    #Generic.solve_fflu(M,bbad)
                    @test_throws DomainError Generic.solve_fflu(M,bbad)
                end
            end
        end
    end

    @testset "Difficult HNF..." begin
        At,xt,bt = ( matrix(ZZ,[-12  19   12  -19   -7;  17   9  -19  -19  -17]),
		   matrix(ZZ, [96  79   81; 18   4  -11; 0 0 0; 0 0 0; 0 0 0]),
		   matrix(ZZ, [135 563 99; 471  -630  -514]))

        
    end
    @testset "Generic.Mat.solve_rational..." begin
        S = ResidueRing(ZZ, 20011*10007)

        for dim = 1:5
            R = MatrixSpace(S, dim, dim)
            U = MatrixSpace(S, dim, rand(1:5))

            M = randmat_with_rank(R, dim, -100:100)
            b = rand(U, -100:100)

            do_test = false
            try
                x, d = solve_rational(M, b)
                do_test = true
            catch e
                if !(e isa ErrorException)
                    rethrow(e)
                end
            end

            if do_test
                @test M*x == d*b
            end
        end

        S, z = PolynomialRing(ZZ, "z")

        for dim = 1:5
            R = MatrixSpace(S, dim, dim)
            U = MatrixSpace(S, dim, rand(1:5))

            M = randmat_with_rank(R, dim, 0:3, -20:20)
            b = rand(U, 0:3, -20:20);

            # Integer division error sometimes, but not always observed.
            x, d = solve_rational(M, b)

            @test M*x == d*b
        end

        R, x = PolynomialRing(QQ, "x")
        K, a = NumberField(x^3 + 3x + 1, "a")

        for dim = 1:5
            S = MatrixSpace(K, dim, dim)
            U = MatrixSpace(K, dim, rand(1:5))

            M = randmat_with_rank(S, dim, 0:2, -100:100)
            b = rand(U, 0:2, -100:100)

            x = solve(M, b)

            @test M*x == b
        end

        R, x = PolynomialRing(ZZ, "x")
        S, y = PolynomialRing(R, "y")

        for dim = 1:5
            T = MatrixSpace(S, dim, dim)
            U = MatrixSpace(S, dim, rand(1:5))

            M = randmat_with_rank(T, dim, 0:2, 0:2, -20:20)
            b = rand(U, 0:2, 0:2, -20:20)

            x, d = solve_rational(M, b)

            @test M*x == d*b
        end

        R, t = PolynomialRing(AbstractAlgebra.JuliaQQ, "t")
        K, a = NumberField(t^3 + 3t + 1, "a")
        S, y = PolynomialRing(K, "y")
        T = MatrixSpace(S, 3, 3)
        U = MatrixSpace(S, 3, 1)

        M = T([3y*a^2 + (y + 1)*a + 2y (5y+1)*a^2 + 2a + y - 1 a^2 + (-a) + 2y; (y + 1)*a^2 + 2y - 4 3y*a^2 + (2y - 1)*a + y (4y - 1)*a^2 + (y - 1)*a + 5; 2a + y + 1 (2y + 2)*a^2 + 3y*a + 3y a^2 + (-y-1)*a + (-y - 3)])
        b = U(permutedims([4y*a^2 + 4y*a + 2y + 1 5y*a^2 + (2y + 1)*a + 6y + 1 (y + 1)*a^2 + 3y*a + 2y + 4], [2, 1]))

        x, d = solve_rational(M, b)

        @test M*x == d*b
    end


    @testset "Generic.Mat.solve_left..." begin
        for R in [ZZ, QQ]
            for iter = 1:4
                for dim = 1:5
                    r = rand(1:5)
                    n = rand(1:5)
                    c = rand(1:5)

                    S = MatrixSpace(R, r, n)
                    U = MatrixSpace(R, n, c)

                    X1 = rand(S, -20:20)
                    M = rand(U, -20:20)

                    B = X1*M
                    X = solve_left(M, X1*M)

                    @test X*M == B
                end
            end
        end

        #=
        R, x = PolynomialRing(QQ, "x")

        for iter = 1:4
            for dim = 1:5
                r = rand(1:5)
                n = rand(1:5)
                c = rand(1:5)

                S = MatrixSpace(R, r, n)
                U = MatrixSpace(R, n, c)

                X1 = rand(S, 1:2, -10:10)
                M = rand(U, 1:2, -10:10)

                B = X1*M
                X = solve_left(M, X1*M)

                @test X*M == B
            end
        end
        =#
    end
#=
    @testset "Generic.Mat.solve_triu..." begin
        R, x = PolynomialRing(QQ, "x")
        K, a = NumberField(x^3 + 3x + 1, "a")

        for dim = 1:10
            S = MatrixSpace(K, dim, dim)
            U = MatrixSpace(K, dim, rand(1:5))

            M = randmat_triu(S, 0:2, -100:100)
            b = rand(U, 0:2, -100:100)

            x = solve_triu(M, b, false)

            @test M*x == b
        end
    end

    @testset "Generic.Mat.solve_left_reduced_triu..." begin
        for iter = 1:40
            n = rand(1:6)
            m = rand(1:n)
            S = MatrixSpace(ZZ, m, n)
            U = MatrixSpace(ZZ, 1, n)

            M = randmat_with_rank(S, rand(1:m), -20:20)
            r = rand(U, -20:20)

            M = hnf(M)

            flag, x = can_solve_left_reduced_triu(r, M)

            @test flag == false || x*M == r
        end
    end
=#
#=
@testset "Generic.Mat.rref..." begin
    S = ResidueRing(ZZ, 20011*10007)
    R = MatrixSpace(S, 5, 5)

    for i = 0:5
        M = randmat_with_rank(R, i, -100:100)

        do_test = false
        r = 0
        d = S(0)
        A = M

        try
            r, d, A = rref(M)
            do_test = true
        catch e
            if !(e isa ErrorException)
                rethrow(e)
            end
        end

        if do_test
            @test r == i
            @test isrref(A)
        end
    end

    S, z = PolynomialRing(ZZ, "z")
    R = MatrixSpace(S, 5, 5)

    for i = 0:5
        M = randmat_with_rank(R, i, 0:3, -20:20)

        r, d, A = rref(M)

        @test r == i
        @test isrref(A)
    end

    R, x = PolynomialRing(QQ, "x")
    K, a = NumberField(x^3 + 3x + 1, "a")
    S = MatrixSpace(K, 5, 5)

    for i = 0:5
        M = randmat_with_rank(S, i, 0:2, -100:100)

        r, A = rref(M)

        @test r == i
        @test isrref(A)
    end

    R, x = PolynomialRing(ZZ, "x")
    S, y = PolynomialRing(R, "y")
    T = MatrixSpace(S, 5, 5)

    for i = 0:5
        M = randmat_with_rank(T, i, 0:2, 0:2, -20:20)

        r, d, A = rref(M)

        @test r == i
        @test isrref(A)
    end
end
=#
end

