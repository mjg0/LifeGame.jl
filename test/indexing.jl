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

        # Ensure that both the underlying value and the index are now correct
        @test lg.grid[I, J] == value
        @test lg[i, j] == true
    end
end
