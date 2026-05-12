"Make sure that correct and pattern are identical"
function checkpattern(correct, pattern)
    identical = all(correct .== pattern)
    if !identical
        printlifediff(correct, correct, pattern; leftlabel="Correct")
    end
    @test identical
end



@testset failfast = true "LifePattern" begin
    rng = MersenneTwister(1)

    for _ in 1:20
        # Create a LifeGrid, a pattern grid, and a LifePattern
        lg = randlifegrid(rng; minsize=(4, 4))
        pattern = rand(rng, Bool, rand(rng, axes(lg, 1)), rand(rng, axes(lg, 2)))
        lp = LifePattern(pattern)

        # Check that setindex works
        I = rand(rng, CartesianIndices(pattern))
        pattern[I] = !pattern[I]
        lp[I] = !lp[I]
        checkpattern(pattern, lp)

        # Check that insertion works
        m, n = size(pattern)
        i = rand(rng, 1:(lastindex(lg, 1)-m+1))
        j = rand(rng, 1:(lastindex(lg, 2)-n+1))
        correct = deepcopy(lg)
        for lpI in CartesianIndices(lp)
            I = lpI + CartesianIndex((i, j)) - oneunit(lpI)
            correct[I] = lp[lpI] || correct[I]
        end
        lgp = insert!(deepcopy(lg), i, j, lp)
        checkpattern(correct, insert!(deepcopy(lg), i, j, pattern))
        checkpattern(correct, lgp)

        # Check that stepping works after insertion
        slow = SlowLifeGrid(lgp; rule=sprint(show, rule(lgp)))
        for _ in 1:10
            step!(lgp)
            step!(slow)
            checkpattern(slow, lgp)
        end
    end
end
