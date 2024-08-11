using Test, LifeGame

include("SlowLifeGrid.jl")



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

    @testset "Indexing of LifeGrid" begin
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

    @testset "Indexing of LifePattern" begin
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

    @testset "step!" begin
        # Test popular rules and rules that used to cause problems for the sparse algorithm
        for rule in ("B3/S23", "B36/S23", "B3678/S34678", "B2/S", "B3/S56", "B37/S357")
            @testset "rule $rule" begin
                for (x, y) in ((1, 1), (4, 5), (15, 61), (35, 63), (326, 256))
                    # Initialize slow and fast grids randomly, with about 1/4 of cells alive
                    grid = rand((false, false, false, true), x, y)
                    slowgrid = SlowLifeGrid(grid; rule=rule)
                    grid_dense_serial    = LifeGrid(grid; rule=rule)
                    grid_sparse_serial   = LifeGrid(grid; rule=rule)
                    grid_dense_parallel  = LifeGrid(grid; rule=rule)
                    grid_sparse_parallel = LifeGrid(grid; rule=rule)
                    # Make sure results are identical over 10 steps
                    for _ in 1:10
                        step!(slowgrid)
                        step!(grid_dense_serial,    dense=true)
                        step!(grid_sparse_serial,   dense=false)
                        step!(grid_dense_parallel,  dense=true,  chunklength=8, parallel=true)
                        step!(grid_sparse_parallel, dense=false, chunklength=8, parallel=true)
                        @test all(slowgrid .== grid_dense_serial
                                           .== grid_sparse_serial
                                           .== grid_dense_parallel
                                           .== grid_sparse_parallel)
                    end
                end
            end
        end
    end

    # @testset "step!" for (x, y) in ((1, 1), (4, 5), (15, 61), (35, 63), (326, 251))
    #     # Test 25 rules, including the default B3/S23
    #     birthrules =    [(3,  ); [(i for i in 1:8 if rand(Bool)) for _ in 1:4]]
    #     survivalrules = [(2, 3); [(i for i in 1:8 if rand(Bool)) for _ in 1:4]]
    #     for birth in birthrules, survival in survivalrules
    #         rule = repr(LifeGame.LifeRule(birth, survival))
    #         # Initialize slow and fast grids randomly, with 1/4 of cells starting alive
    #         grid = rand((false, false, false, true), x, y)
    #         slowgrid = SlowLifeGrid(grid; rule=rule)
    #         fastgrid_dense_serial    = LifeGrid(grid; rule=rule)
    #         fastgrid_sparse_serial   = LifeGrid(grid; rule=rule)
    #         fastgrid_dense_parallel  = LifeGrid(grid; rule=rule)
    #         fastgrid_sparse_parallel = LifeGrid(grid; rule=rule)
    #         # Make sure results are identical over 5 steps
    #         for _ in 1:5
    #             step!(slowgrid)
    #             step!(fastgrid_dense_serial,    sparse=false)
    #             step!(fastgrid_sparse_serial,   sparse=true)
    #             step!(fastgrid_dense_parallel,  sparse=false, chunklength=8, parallel=true)
    #             step!(fastgrid_sparse_parallel, sparse=true,  chunklength=8, parallel=true)
    #             @test all(slowgrid .== fastgrid_dense_serial
    #                                .== fastgrid_sparse_serial
    #                                .== fastgrid_dense_parallel
    #                                .== fastgrid_sparse_parallel)
    #         end
    #     end
    # end
end # @testset
