using LifeGame # LifeGame.LifeRule, LifeGame.step!



# Basic life grid implementation against which to test LifeGrid
mutable struct SlowLifeGrid <: AbstractMatrix{Bool}
    grid::Matrix{Bool} # the grid itself
    next::Matrix{Bool} # area to store intermediate results
    birthsums::Vector{Int}    # list of neighbor sums that lead to cell birth
    survivalsums::Vector{Int} # list of neighbor sums that allow cell survival

    function SlowLifeGrid(height, width; rule="B3/S23")
        birthsums, survivalsums = LifeGame.rulesums(LifeGame.LifeRule(rule))
        return new(zeros(height, width), zeros(height, width), birthsums, survivalsums)
    end

    SlowLifeGrid(grid; kw...) = SlowLifeGrid(size(grid)...; kw...) .= grid
end



# Implement AbstractArray interface for SlowLifeGrid
Base.size(lg::SlowLifeGrid) = size(lg.grid)

Base.@propagate_inbounds function Base.getindex( lg::SlowLifeGrid, x...)
    return getindex( lg.grid, x...)
end

Base.@propagate_inbounds function Base.setindex!(lg::SlowLifeGrid, x...)
    return setindex!(lg.grid, x...)
end



# Update a SlowLifeGrid
function LifeGame.step!(lg::SlowLifeGrid)
    # Iterate over the whole grid
    region = CartesianIndices(lg.grid)
    @inbounds @simd for I in region
        # Sum the living neighbors of the current cell
        neighborhood = max(first(region), I-oneunit(I)):min(last(region), I+oneunit(I))
        neighborsum = sum(@view lg.grid[neighborhood]) - lg.grid[I]

        survival = neighborsum in lg.survivalsums
        birth    = neighborsum in lg.birthsums

        lg.next[I] = lg.grid[I] && survival || birth
    end

    # Swap current and next grids and return the updated SlowLifeGrid
    lg.grid, lg.next = lg.next, lg.grid
    return lg
end
