@testset "LifeGrid constructor" begin
    # Basics
    @test (
        LifeGrid(2, 3) ==
        LifeGrid(2, 3; rule="B3/S23") ==
        LifeGrid([0 0 0; 0 0 0]) ==
        LifeGrid(zeros(Bool, 2, 3)) ==
        LifeGrid(BitArray([0 0 0; 0 0 0]))
    )

    # Indexing from newly constructed grid
    @test all(
        LifeGrid([0 1 1 0; 1 1 0 0; 0 1 0 1]) .==
        [false true true false; true true false false; false true false true],
    )

    # White box test of a known good case
    grid = [
        0 1 0 1 0 1 0
        0 0 1 1 0 0 1
        0 1 1 0 1 1 0
        1 0 0 0 1 0 1
        0 1 0 0 1 1 0
    ]
    lg = LifeGrid(grid, CType=UInt8, HType=UInt8)
    @test all(lg .== grid)
    @test all(
        lg.grid[1:5, :] .== [
            0b00101010 0b00000000
            0b00011000 0b01000000
            0b00110110 0b00000000
            0b01000100 0b01000000
            0b00100110 0b00000000
        ],
    )
    @test all(lg.halos.currentleft[:] .== [0b00010000, 0b01010000])
    @test all(lg.halos.currentright[:] .== [0b10101000, 0b00000000])
end
