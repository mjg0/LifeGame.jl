using LifeGameMPI, LifeGame, MPI

MPI.Init()

# Create a random grid with a trillion cells
lg =  LifeGridMPI(250_000, 250_000)
lg.lifegrid.grid[begin+1:end-1,begin+1:end-1] .= rand(eltype(lg.lifegrid.grid), size(lg.lifegrid.grid).-2)

# Compile and sync up
LifeGameMPI.step!(lg) # compile first
MPI.Barrier(MPI.COMM_WORLD)

# See how long it takes to step 100 times
@time for _ in 1:25
    LifeGameMPI.step!(lg)
end
