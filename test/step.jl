@testset failfast = true "step!" begin
    rng = MersenneTwister(1)

    # First run a small test that will print a helpful message on failure
    lgsmall = LifeGrid(rand(rng, Bool, 15, 20); CType=UInt8, HType=UInt8)
    slowsmall = SlowLifeGrid(lgsmall)
    for _ in 1:10
        prev = deepcopy(lgsmall)
        step!(lgsmall)
        step!(slowsmall)
        identical = all(lgsmall .== slowsmall)
        if !identical
            printlifediff(prev, slowsmall, lgsmall; leftlabel="Previous")
        end
        @test identical
    end

    # Run many tests with different sizes and rules
    for _ in 1:100
        lg = randlifegrid(rng)
        threadmode = rand(rng, (serial, parallel))
        slow = SlowLifeGrid(lg)
        for _ in 1:10 # stepping costs almost nothing compared to all the compiling
            @test all(step!(lg, threadmode) .== step!(slow))
        end
    end
end