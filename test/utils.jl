"Return a random life rule string (e.g. 'B3/S23')"
function randrule(rng=default_rng())
    randrulevec() = [i for i = 1:8 if rand(rng, Bool)]

    R = LifeGame.LifeRule(randrulevec(), randrulevec())

    return sprint(show, R)
end



"Return a random LifeGrid"
function randlifegrid(rng=default_rng(); minsize=(1, 1), maxsize=(300, 300))
    R = randrule(rng)

    C, H = (rand(rng, (UInt8, UInt16, UInt32, UInt64)) for _ in 1:2)

    lgsize = (rand(rng, a:b) for (a, b) in zip(minsize, maxsize))

    return LifeGrid(lgsize...; rule=R, CType=C, HType=H)
end



"""
    printlinediff(left, correct, actual; leftlabel)

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