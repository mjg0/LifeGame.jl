function checkpattern(correct, pattern)
    identical = all(correct .== pattern)
    if !identical && max(size(correct)...) < 30
        printlifediff(correct, correct, pattern; leftlabel="Correct")
    end
    @test identical
end



@testset failfast = true "LifePattern" begin
    rng = MersenneTwister(1)

    for _ in 1:20
        # Create a LifeGrid, a pattern grid, and a LifePattern
        lg = randlifegrid(rng; minsize=(4, 4), maxsize=(25, 25))
        pattern = rand(rng, Bool, rand(rng, axes(lg, 1)), rand(rng, axes(lg, 2)))
        lp = LifePattern(pattern)

        # Check that setindex works
        I = rand(rng, CartesianIndices(pattern))
        pattern[I] = !pattern[I]
        lp[I] = !lp[I]
        checkpattern(pattern, lp)

        # Check that insertion works
        m, n = size(pattern)
        i = rand(rng, 1:(lastindex(lg, 1)-m))
        j = rand(rng, 1:(lastindex(lg, 2)-n))
        correct = deepcopy(lg)
        for lpI in CartesianIndices(lp)
            I = lpI + CartesianIndex((i, j)) - oneunit(lpI)
            correct[I] = lp[lpI] || correct[I]
        end
        lgp = insert!(deepcopy(lg), i, j, lp)
        checkpattern(correct, insert!(deepcopy(lg), i, j, pattern))
        checkpattern(correct, lgp)

        # Check that stepping works after insertion
        slow = SlowLifeGrid(lgp)
        for _ in 1:10
            checkpattern(slow, lgp)
            step!(lgp)
            step!(slow)
            #@test all(step!(lgp) .== step!(slow))
        end
    end
end
