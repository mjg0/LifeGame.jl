using Test, Crayons, Random, LifeGame

include("SlowLifeGrid.jl")



# Print what went wrong when a LifeGrid was updated incorrectly
function printlifediff(left, correct, actual; leftlabel)
    # Basic information
    m, n = size(correct)
    x, y = Tuple(findfirst(correct .!= actual))
    println("Failure with $m×$n grid, first difference at cell ($x, $y)")

    # Short-circuit if the grid is too big
    if max(size(correct)...) > 20
        return
    end

    # Header
    println("Incorrect result (bad cells in red):")
    printedwidth = 2*size(left, 2)+4
    print(leftlabel, repeat(' ', max(0, printedwidth-length(leftlabel))))
    println("Actual", repeat(' ', max(0, printedwidth-length("Actual"))))

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
                print(Crayon(foreground = :red), cell, Crayon(reset = true))
            else
                print(cell)
            end
        end

        println()
    end
end



# step! LifeGrid a few times 
function teststep!(rule, size; CType = nothing)
    m, n = size
    @testset failfast=true "rule $rule, size $m×$n" begin
        for parallel in (false, true)
            # Make results reproducible
            rng = MersenneTwister(1)

            # Initialize
            grid = rand(rng, (true, false, false, false), m, n)
            slow = SlowLifeGrid(grid; rule = rule)
            lg = if CType === nothing
                LifeGrid(grid; rule = rule)
            else
                LifeGrid(grid; rule = rule, CType = CType)
            end

            # Function to step then test
            function stepandtest()
                prev = deepcopy(slow)
                step!(slow)
                step!(lg; parallel = parallel)
                gridequal = all(slow .== lg)
                if !gridequal
                    println("Failed with rule $rule")
                    printlifediff(prev, slow, lg; leftlabel = "Previous")
                end
                @test gridequal
            end

            # Test once, flip a random cell, then test again a few more times
            stepandtest()
            for _ = 1:3
                I = rand(rng, CartesianIndices(lg))
                slow[I] = !slow[I]
                lg[I] = !lg[I]
                stepandtest()
            end
        end
    end
end



# Tests
@testset "LifeGame" begin
    @testset "LifeGrid construction" begin
        @test (
            LifeGrid(2, 3) ==
            LifeGrid(2, 3; rule = "B3/S23") ==
            LifeGrid([0 0 0; 0 0 0]) ==
            LifeGrid(zeros(Bool, 2, 3)) ==
            LifeGrid(BitArray([0 0 0; 0 0 0]))
        )
        @test all(
            LifeGrid([0 1 1 0; 1 1 0 0; 0 1 0 1]) .==
            [false true true false; true true false false; false true false true],
        )
    end

    @testset "LifeGrid indexing" begin
        lg = LifeGrid(5, 200)
        # Check indexing in several cases, especially at borders between clusters
        for (i, j, I, J, value) in (
            (1, 1, 1, 1, one(UInt64) << 62),
            (2, 2, 2, 1, one(UInt64) << 61),
            (3, 62, 3, 1, one(UInt64) << 1),
            (4, 63, 4, 2, one(UInt64) << 62),
            (5, 125, 5, 3, one(UInt64) << 62),
            (5, 126, 5, 3, UInt64(0x03) << 61),
        )
            # The grid starts zeroed, so the cell should start false
            @test lg[i, j] == false
            # Both numbers and booleans should be accepted by setindex!
            lg[i, j] = 1
            lg[i, j] = true
            # Ensure that both the underlying value and the index are correct
            @test lg.grid[I+1, J+1] == value
            @test lg[i, j] == true
        end
    end

    @testset failfast=true "LifePattern" begin
        rng = MersenneTwister(1)
        sizes = (5, 20, 80, 150)
        CTypes = (UInt8, UInt16, UInt32, UInt64)

        for M in sizes, N in sizes
            lg = LifeGrid(rand(rng, Bool, 500, 500), CType = rand(rng, CTypes))

            # Construction
            pattern = rand(rng, Bool, M, N)
            lp = LifePattern(pattern)

            # Indexing
            I = rand(rng, CartesianIndices(pattern))
            pattern[I] = !pattern[I]
            lp[I] = !lp[I]
            identical = all(lp .== pattern)
            if !identical
                printlifediff(pattern, pattern, lp; leftlabel = "Correct")
            end
            @test identical

            # Insertion
            i = rand(rng, 1:(lastindex(lg, 1)-M))
            j = rand(rng, 1:(lastindex(lg, 2)-N))
            correct = Array(lg)
            for lpI in CartesianIndices(lp)
                I = lpI+CartesianIndex((i, j))-oneunit(lpI)
                correct[I] = lp[lpI] || correct[I]
            end
            lgwithp1 = insert!(deepcopy(lg), i, j, pattern)
            lgwithp2 = insert!(deepcopy(lg), i, j, lp)
            identical = all(correct .== lgwithp1 .== lgwithp2)
            if !identical
                printlifediff(correct, correct, lgwithp2; leftlabel = "Correct")
            end
            @test identical

            # Check halos by stepping
            stepped = step!(LifeGrid(correct))
            step!(lgwithp1)
            step!(lgwithp2)
            identical = all(stepped .== lgwithp1 .== lgwithp2)
            if !identical
                printlifediff(correct, stepped, lgwithp2; leftlabel = "Previous")
            end
        end
    end



    @testset "updatedcluster" begin
        # Test with some fixed values that have been calculated by hand
        rule = LifeGame.LifeRule("B3/S23")
        for (above, middle, below, result) in (
            (0b1100, 0b1000, 0b0000, 0b1100),
            (0b0100, 0b0100, 0b0100, 0b1110),
            (0b0010, 0b1010, 0b0110, 0b0011),
            (0b1000, 0b0110, 0b1100, 0b0010),
        )
            @test LifeGame.updatedcluster(above, middle, below, rule) == result
        end
    end

    @testset "updatedchunkhalos" begin
        chunk = [
            0b0100000000000000,
            0b0000000000000010,
            0b0100000000000010,
            0b0000000000000000,
            0b0100000000000010,
            0b0000000000000010,
            0b0100000000000000,
            0b0000000000000010,
        ]
        lhalo = 0b10101010
        rhalo = 0b01101101

        lactual, ractual = LifeGame.updatedchunkhalos(chunk, UInt8)

        @test lhalo == lactual
        @test rhalo == ractual
    end

    @testset "updatehalos!" begin
        lhalo = 0b10101010
        rhalo = 0b01101101
        in = zeros(UInt16, 8)
        out = zeros(UInt16, 10)
        correct = [
            0b0000000000000000,
            0b1000000000000000,
            0b0000000000000001,
            0b1000000000000001,
            0b0000000000000000,
            0b1000000000000001,
            0b0000000000000001,
            0b1000000000000000,
            0b0000000000000001,
            0b0000000000000000,
        ]

        LifeGame.updatehalos!(out, in, lhalo, rhalo)

        @test all(out .== correct)
    end

    @testset "LifeGrid construction" begin
        grid = [
            0 1 0 1 0 1 0
            0 0 1 1 0 0 1
            0 1 1 0 1 1 0
            1 0 0 0 1 0 1
            0 1 0 0 1 1 0
        ]
        lg = LifeGrid(grid, CType = UInt8, HType = UInt8)
        @test all(lg .== grid)
        @test all(
            lg.grid[2:6, 2:3] .== [
                0b00101010 0b00000000
                0b00011000 0b01000000
                0b00110110 0b00000000
                0b01000100 0b01000000
                0b00100110 0b00000000
            ],
        )
        @test all(lg.halos.currentleft[2:3] .== [0b00010000, 0b01010000])
        @test all(lg.halos.currentright[2:3] .== [0b10101000, 0b00000000])
    end

    @testset "step!" begin
        rng = MersenneTwister(1)
        CTypes = (UInt8, UInt16, UInt32, UInt64)
        rules = ("B3/S23", "B2/S", "B35678/S5678", "B36/S23", "B234/S", "B345/S5")

        # Numerous tests of grids of size close to the edge of clusters and chunks
        for height in (3, 7, 8, 9, 40, 64, 65, 128, 129, 150, 300)
            for width in (6, 7, 14, 15, 45, 62, 63, 90, 124, 125, 200, 400)
                teststep!(rand(rng, rules), (height, width); CType = rand(rng, CTypes))
            end
        end
    end
end # @testset
