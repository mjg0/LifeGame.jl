using Test, LifeGame



# Basic life grid implementation against which to test LifeGrid
mutable struct SlowLifeGrid <: AbstractMatrix{Bool}
    grid::Matrix{Bool}
    next::Matrix{Bool}

    SlowLifeGrid(height, width) = new(zeros(height, width), zeros(height, width))

    SlowLifeGrid(grid) = SlowLifeGrid(size(grid)...) .= grid
end



# Implement AbstractArray interface for SlowLifeGrid
Base.size(lg::SlowLifeGrid) = size(lg.grid)

Base.@propagate_inbounds Base.getindex( lg::SlowLifeGrid, x...) = getindex( lg.grid, x...)

Base.@propagate_inbounds Base.setindex!(lg::SlowLifeGrid, x...) = setindex!(lg.grid, x...)



# Update a SlowLifeGrid
function LifeGame.step!(lg::SlowLifeGrid)
    R = CartesianIndices(lg.grid)
    @inbounds @simd for I in R
        # Sum the living neighbors (including the cell in question) of the current cell
        alivecount = 0
        for cell in max(first(R), I-oneunit(I)):min(last(R), I+oneunit(I))
            if lg.grid[cell]
                alivecount += 1
            end
        end
        # Update the cell
        lg.next[I] = alivecount == 3 || alivecount == 4 && lg.grid[I]
    end

    # Swap current and next grids and return the updated SlowLifeGrid
    lg.grid, lg.next = lg.next, lg.grid
    return lg
end



# Tests
@testset "LifeGame" begin
    @testset "LifeGrid construction" begin
        @test (LifeGrid(2, 3) == LifeGrid([0 0 0; 0 0 0])
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
        for (above, middle, below, result) in ((0b1100, 0b1000, 0b0000, 0b1100),
                                               (0b0100, 0b0100, 0b0100, 0b1110),
                                               (0b0010, 0b1010, 0b0110, 0b0011),
                                               (0b1000, 0b0110, 0b1100, 0b0010))
            @test LifeGame.updatedcluster(above, middle, below) == result
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

    @testset "step!" for (x, y) in ((1, 1), (4, 5), (15, 61), (35, 63), (326, 251))
        # Initialize slow and fast grids randomly, with 1/4 of cells starting alive
        grid = rand((false, false, false, true), x, y)
        slowgrid = SlowLifeGrid(grid)
        fastgrid_dense_serial    = LifeGrid(grid)
        fastgrid_sparse_serial   = LifeGrid(grid)
        fastgrid_dense_parallel  = LifeGrid(grid)
        fastgrid_sparse_parallel = LifeGrid(grid)
        # Make sure results are identical over 10 steps
        for _ in 1:10
            step!(slowgrid)
            step!(fastgrid_dense_serial,    sparse=false)
            step!(fastgrid_sparse_serial,   sparse=true)
            step!(fastgrid_dense_parallel,  sparse=false, chunklength=8, parallel=true)
            step!(fastgrid_sparse_parallel, sparse=true,  chunklength=8, parallel=true)
            @test all(slowgrid .== fastgrid_dense_serial
                               .== fastgrid_sparse_serial
                               .== fastgrid_dense_parallel
                               .== fastgrid_sparse_parallel)
        end
    end
end # @testset