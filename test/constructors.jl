@testset "LifeRule constructor" begin
    rule = LifeGame.LifeRule("B024/S056")
    @test repr(rule) == "B024/S056"
    @test LifeGame.rulesums(rule) == ([0, 2, 4], [0, 5, 6])
    @test LifeGame.Rule(0, 8) == 0x0101
    @test_throws ArgumentError LifeGame.Rule(9)
end



@testset "LifeGrid constructor" begin
    rng = MersenneTwister(1)

    # Get C and H from an exising LifeGrid
    ctype(::LifeGrid{R,C,H,T,W}) where {R,C,H,T,W} = C
    htype(::LifeGrid{R,C,H,T,W}) where {R,C,H,T,W} = H

    # Default C and H type for a given grid size
    defaultctype(lg::LifeGrid) = LifeGame.smallestuint(size(lg, 2) + 2)
    defaulthtype(lg::LifeGrid) = LifeGame.smallestuint(size(lg, 1))

    @testset "constructor" begin
        # Default rule is B3/S23
        @test repr(rule(LifeGrid(3, 4))) == "B3/S23"

        for _ in 1:20
            lg = randlifegrid(rng)

            # Identical excepting HType
            lg1 = LifeGrid(lg; rule=rule(lg), CType=ctype(lg))
            lg2 = LifeGrid(size(lg)...; rule=rule(lg), CType=ctype(lg))
            @test sum(lg2) == 0
            lg2 .= lg1
            @test all(lg1 .== lg2 .== lg)
            @test htype(lg1) == htype(lg2) == defaulthtype(lg)
            @test ctype(lg1) == ctype(lg2) == ctype(lg1)
            @test rule(lg1) == rule(lg2) == rule(lg)

            # Identical excepting CType
            lg3 = LifeGrid(lg; rule=rule(lg), HType=htype(lg))
            lg4 = LifeGrid(size(lg)...; rule=rule(lg), HType=htype(lg))
            @test sum(lg4) == 0
            lg4 .= lg3
            @test all(lg3 .== lg4 .== lg)
            @test htype(lg3) == htype(lg4) == htype(lg)
            @test ctype(lg3) == ctype(lg4) == defaultctype(lg)
            @test rule(lg3) == rule(lg4) == rule(lg)

            # CType and HType both unspecified
            lg5 = LifeGrid(lg; rule=rule(lg))
            lg6 = LifeGrid(size(lg)...; rule=rule(lg))
            @test sum(lg6) == 0
            lg6 .= lg5
            @test all(lg5 .== lg6 .== lg)
            @test htype(lg5) == htype(lg6) == defaulthtype(lg)
            @test ctype(lg5) == ctype(lg6) == defaultctype(lg)
            @test rule(lg5) == rule(lg6) == rule(lg)
        end
    end
end



@testset "LifePattern constructor" begin
    rng = MersenneTwister(1)

    for _ in 1:20
        # Build two patterns, an empty one and a random one
        grid = rand(rng, Bool, rand(rng, 1:200), rand(rng, 1:200))
        pattern = LifePattern(grid)
        emptypattern = LifePattern(size(pattern)...)

        # Sizes
        @test all(size(emptypattern) .== size(pattern) .== size(grid))

        # Values
        @test sum(emptypattern) == 0
        @test all(pattern .== grid)

        # Indexing
        I = rand(rng, CartesianIndices(grid))
        grid[I] = !grid[I]
        pattern[I] = !pattern[I]
        @test all(pattern .== grid)
    end
end
