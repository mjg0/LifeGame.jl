@testset "LifeGrid I/O" begin
    rng = MersenneTwister(1)

    for _ = 1:20
        # Create a random LifeGrid
        lg = randlifegrid(rng)

        # Write it to a stream
        io = IOBuffer()
        write(io, lg)

        # Read in a copy
        seekstart(io)
        lgc = LifeGrid(io)

        # Ensure they're identical now...
        @test all(lg .== lgc)

        # ...and after stepping
        for _ in 1:10
            @test all(step!(lg) .== step!(lgc))
        end
    end
end