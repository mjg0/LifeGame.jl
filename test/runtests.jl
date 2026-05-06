using Test, Crayons, Random, LifeGame

include("SlowLifeGrid.jl")



function printlifediff(previous, correct, actual)
    # Header
    m, n = size(actual)
    println("Incorrect result for $m×$n $(rule(actual)) grid")
    println("Bad cells are printed in red")
    printedwidth = 2*size(previous, 2)+4
    print(  "Previous", repeat(' ', max(0, printedwidth-8)))
    println("Computed", repeat(' ', max(0, printedwidth-8)))

    # Representation of whether the cell is on or off
    rep(onoroff) = onoroff ? "* " : "- "

    for i in axes(correct, 1)
        # Print the previous grid
        for j in axes(correct, 2)
            print(rep(previous[i,j]))
        end
        print("    ")

        # Print the actual grid, highlighting incorrect cells
        for j in axes(correct, 2)
            cell = rep(actual[i,j])
            if actual[i,j] != correct[i,j]
                print(Crayon(foreground=:red), cell, Crayon(reset=true))
            else
                print(cell)
            end
        end

        println()
    end
end



function teststep!(rule, size; verbose=false, CType=nothing)
    m, n = size
    @testset failfast=true "rule $rule, size $m×$n" begin
        for parallel in (false, true)
            # Make results reproducible
            rng = MersenneTwister(1)

            # Initialize
            grid = rand(rng, (true, false, false, false), m, n)
            slow = SlowLifeGrid(grid; rule=rule)
            lg = if CType === nothing
                LifeGrid(grid; rule=rule)
            else
                LifeGrid(grid; rule=rule, CType=CType)
            end

            # Function to step then test
            function stepandtest()
                prev = deepcopy(slow)
                step!(slow)
                step!(lg; parallel=parallel)
                gridequal = all(slow .== lg)
                if !gridequal
                    println("Failed for rule $rule, size $m×$n")
                    println("First different index: $(findfirst(==(false), slow .== lg))")
                    if verbose
                        printlifediff(prev, slow, lg)
                    end
                end
                @test gridequal
            end

            # Test once, flip a random cell, then test again a few more times
            stepandtest()
            for _ in 1:3
                I = rand(rng, CartesianIndices(lg))
                slow[I] = !slow[I]
                lg[  I] = !lg[  I]
                stepandtest()
            end
        end
    end
end



# Tests
@testset "LifeGame" begin
    @testset "LifeGrid construction" begin
        @test (LifeGrid(2, 3) == LifeGrid(2, 3; rule="B3/S23")
                              == LifeGrid([0 0 0; 0 0 0])
                              == LifeGrid(zeros(Bool, 2, 3))
                              == LifeGrid(BitArray([0 0 0; 0 0 0])))
        @test all(LifeGrid([0 1 1 0; 1 1 0 0; 0 1 0 1]) .==
                  [false true true false; true true false false; false true false true])
    end

    @testset "LifeGrid indexing" begin
        lg = LifeGrid(5, 200)
        # Check indexing in several cases, espcially at borders between clusters
        for (i, j, I, J, value) in ((1, 1,   1, 1, one(UInt64)  << 62),
                                    (2, 2,   2, 1, one(UInt64)  << 61),
                                    (3, 62,  3, 1, one(UInt64)  << 1 ),
                                    (4, 63,  4, 2, one(UInt64)  << 62),
                                    (5, 125, 5, 3, one(UInt64)  << 62),
                                    (5, 126, 5, 3, UInt64(0x03) << 61))
            # The grid starts zeroed, so the cell should start false
            @test lg[i,j] == false
            # Both numbers and booleans should be accepted by setindex!
            lg[i,j] = 1
            lg[i,j] = true
            # Ensure that both the underlying value and the index are correct
            @test lg.grid[I+1,J+1] == value
            @test lg[i,j] == true
        end
    end

    @testset "LifePattern indexing" begin
        pat = LifePattern([0 1 0
                           0 0 1
                           1 1 1])
        # Test a few places, especially the beginning and end of the pattern
        @test pat[1,1] == false
        @test pat[1,2] == true
        @test pat[3,2] == true
        # Flip a few values, both with a bool and with a number
        pat[3,3] = false
        pat[2,1] = 1
        # Test that the swaps are reflected by getindex
        @test pat[3,3] == false
        @test pat[2,3] == true
    end

    @testset "updatedcluster" begin
        # Test with some fixed values that have been calculated by hand
        rule = LifeGame.LifeRule("B3/S23")
        for (above, middle, below, result) in ((0b1100, 0b1000, 0b0000, 0b1100),
                                               (0b0100, 0b0100, 0b0100, 0b1110),
                                               (0b0010, 0b1010, 0b0110, 0b0011),
                                               (0b1000, 0b0110, 0b1100, 0b0010))
            @test LifeGame.updatedcluster(above, middle, below, rule) == result
        end
    end

    @testset "LifePattern insertion" begin
        lg = LifeGrid(10, 100)
        # Insert a few patterns into a grid, especially at cluster borders
        for (pattern, i, j) in (([1 0 1], 1, 1),
                                ([1 1 0 0 1 1 0 0 1
                                  0 0 1 1 1 0 0 1 1], 2, 1),
                                ([1 1 1 0 1
                                  1 0 0 1 0], 4, 60),
                                (rand(Bool, 3, 70), 6, 20))
            insert!(lg, i, j, LifePattern(pattern))
            x, y = size(pattern)
            # The inserted pattern should be present at the correct place in the LifeGrid
            @test all(lg[i:i+x-1,j:j+y-1] .== pattern)
        end
        # Make sure that an exception is thrown when out-of-bounds access is attempted
        @test_throws BoundsError insert!(LifeGrid(4, 5), 2, 4, LifePattern([1 1 1]))
    end

    @testset "updatedchunkhalos" begin
        chunk = [0b0100000000000000,
                 0b0000000000000010,
                 0b0100000000000010,
                 0b0000000000000000,
                 0b0100000000000010,
                 0b0000000000000010,
                 0b0100000000000000,
                 0b0000000000000010]
        lhalo = 0b10101010
        rhalo = 0b01101101

        lactual, ractual = LifeGame.updatedchunkhalos(chunk, UInt8)

        @test lhalo == lactual
        @test rhalo == ractual
    end

    @testset "updatehalos!" begin
        lhalo = 0b10101010
        rhalo = 0b01101101
        in  = zeros(UInt16, 8)
        out = zeros(UInt16, 10)
        correct = [0b0000000000000000,
                   0b1000000000000000,
                   0b0000000000000001,
                   0b1000000000000001,
                   0b0000000000000000,
                   0b1000000000000001,
                   0b0000000000000001,
                   0b1000000000000000,
                   0b0000000000000001,
                   0b0000000000000000]
        
        LifeGame.updatehalos!(out, in, lhalo, rhalo)

        @test all(out .== correct)
    end

    @testset "LifeGrid construction" begin
        grid = [0 1 0 1 0 1   0
                0 0 1 1 0 0   1
                0 1 1 0 1 1   0
                1 0 0 0 1 0   1
                0 1 0 0 1 1   0]
        lg = LifeGrid(grid, CType=UInt8, HType=UInt8)
        @test all(lg .== grid)
        @test all(lg.grid[2:6,2:3] .== [0b00101010 0b00000000
                                        0b00011000 0b01000000
                                        0b00110110 0b00000000
                                        0b01000100 0b01000000
                                        0b00100110 0b00000000])
        @test all(lg.lefthalos[ 1][2:3] .== [0b00010000, 0b01010000])
        @test all(lg.righthalos[1][2:3] .== [0b10101000, 0b00000000])
    end

    @testset "step!" begin
        # A few small grid tests with verbose output on failure
        for rule in ("B3/S23", "B2/S", "B35678/S5678")
            for size in ((3, 4), (8, 6), (9, 7))
                teststep!(rule, size; verbose=true, CType=UInt8)
            end
        end

        # Numerous tests of grids of size close to the edge of clusters and chunks
        for rule in ("B36/S23", "B234/S", "B345/S5")
            for height in (8, 40, 64, 65, 128, 129, 150, 300)
                for width in (6, 45, 62, 63, 90, 124, 125, 200, 400)
                    teststep!(rule, (height, width); verbose=false)
                end
            end
        end
   end
end # @testset
