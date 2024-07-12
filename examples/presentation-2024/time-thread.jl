using LifeGame

# Create a random grid with a trillion cells
lg =  LifeGrid(250_000, 250_000)
lg.grid[begin+1:end-1,begin+1:end-1] .= rand(eltype(lg.grid), size(lg.grid).-2)

# See how long it takes to step 100 times
step!(lg) # compile first
@time for _ in 1:25
    step!(lg)
end
