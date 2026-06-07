"""
    randlifegrid(rng)

Return a randomly initialized `LifeGrid` of random size with a random `LifeRule`.
"""
function randlifegrid(rng)
    # Cluster and halo types
    utypes = (UInt8, UInt16, UInt32, UInt64)
    C = rand(rng, utypes)
    H = rand(rng, utypes)

    # Cluster size, biased toward sizes likely to hit fencepost errors
    gridsize = (rand(rng, (1, 2, 3, n-1, n, n+1, 2n-1, 2n, 2n+1, rand(rng, 3n:4n)))
                for n in (8 * sizeof(H), 8 * sizeof(C) - 2))

    # Rule
    birth, survival = ntuple(_ -> LifeGame.Rule{rand(rng, UInt8)}, 2)
    rule = LifeGame.LifeRule{birth,survival}()

    # Build and randomize the LifeGrid
    lg = LifeGrid(gridsize...; rule=repr(rule), CType=C, HType=H)
    rand!(rng, lg)

    return lg
end




"""
    printlifediff(left, correct, actual; leftlabel)

Print highlighted differences between `correct` and `actual`

`left` is printed on the left side, with `leftlabel` as a header. `actual` is printed on the
right side, with differences between `actual` and `correct` printed in red.
"""
function printlifediff(left, correct, actual; leftlabel)
    # Basic information
    m, n = size(correct)
    x, y = Tuple(findfirst(correct .!= actual))
    println("Failure with $m×$n grid, first difference at cell ($x, $y)")

    # Short-circuit if the grid is too big
    if max(size(correct)...) > 30
        return
    end

    # Header
    println("Incorrect result (bad cells in red):")
    printedwidth = 2 * size(left, 2) + 4
    print(leftlabel, repeat(' ', max(0, printedwidth - length(leftlabel))))
    println("Actual", repeat(' ', max(0, printedwidth - length("Actual"))))

    # Representation of whether the cell is on or off
    rep(onoroff) = onoroff ? "* " : "- "

    for i in axes(correct, 1)
        # Print the left grid
        for j in axes(correct, 2)
            print(rep(left[i, j]))
        end
        print("    ")

        # Print the actual grid, highlighting incorrect cells
        for j in axes(correct, 2)
            cell = rep(actual[i, j])
            if actual[i, j] != correct[i, j]
                print(Crayon(foreground=:red), cell, Crayon(reset=true))
            else
                print(cell)
            end
        end

        println()
    end
end



"""
    stepandcheck!(reference, grid, rng, mutate! = x -> nothing)

`mutate!` `reference` and `grid`, `step!` them, and test that they remain identical.

`mutate!` is a function that takes a grid as its argument and optionally mutates it.
`reference` is the known correct grid, for which `mutate!` and `step!` must both work
correctly. `grid` is the grid to be tested, and `rng` is a random number generator.
"""
function stepandcheck!(reference, grid, rng, mutate! = x -> nothing)
    mutate!(reference)
    before = deepcopy(reference)
    step!(reference)

    mutate!(grid)
    threadmode = rand(rng, (serial, parallel))

    step!(grid, threadmode)

    correct = all(grid .== reference)
    if !correct
        println("Failure with rule $(rule(grid))")
        printlifediff(before, reference, grid; leftlabel="Previous")
    end
    @test correct
end