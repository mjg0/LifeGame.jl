@testset "updatedcluster" begin
    # Test with some fixed values that have been calculated by hand
    R = LifeGame.LifeRule("B3/S23")

    for (above, middle, below, result) in (
        (0b1100, 0b1000, 0b0000, 0b1100),
        (0b0100, 0b0100, 0b0100, 0b1110),
        (0b0010, 0b1010, 0b0110, 0b0011),
        (0b1000, 0b0110, 0b1100, 0b0010),
    )
        @test LifeGame.updatedcluster(above, middle, below, R) == result
    end
end
