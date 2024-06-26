using LifeGame # since we overload LifeGame.step!



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
    # Iterate over the whole grid
    R = CartesianIndices(lg.grid)
    @inbounds @simd for I in R
        # Sum the living neighbors (including the cell in question) of the current cell
        alivecount = zero(UInt8)
        for neighbor in max(first(R), I-oneunit(I)):min(last(R), I+oneunit(I))
            alivecount += reinterpret(UInt8, lg.grid[neighbor])
        end

        # Update the cell
        lg.next[I] = alivecount == 3 || alivecount == 4 && lg.grid[I]
    end

    # Swap current and next grids and return the updated SlowLifeGrid
    lg.grid, lg.next = lg.next, lg.grid
    return lg
end
