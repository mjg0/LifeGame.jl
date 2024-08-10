using LifeGame # LifeGame.LifeRule, LifeGame.step!



# Basic life grid implementation against which to test LifeGrid
mutable struct SlowLifeGrid <: AbstractMatrix{Bool}
    grid::Matrix{UInt8} # the grid itself
    next::Matrix{UInt8} # area to store intermediate results
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
    return reinterpret(Bool, getindex( lg.grid, x...))
end

Base.@propagate_inbounds function Base.setindex!(lg::SlowLifeGrid, x...)
    return reinterpret(Bool, setindex!(lg.grid, x...))
end



# Update a SlowLifeGrid
function LifeGame.step!(lg::SlowLifeGrid)
    # Iterate over the whole grid
    region = CartesianIndices(lg.grid)
    @inbounds @simd for I in region
        # Sum the living neighbors (excluding the cell in question) of the current cell
        neighborhood = max(first(region), I-oneunit(I)):min(last(region), I+oneunit(I))
        neighborsum = sum(n->lg.grid[n], neighborhood) - lg.grid[I]

        # Kill cells that shouldn't survive
        lg.next[I] = neighborsum in lg.survivalsums ? lg.grid[I] : 0x00

        # Spawn cells that should be born
        lg.next[I] = neighborsum in    lg.birthsums ? 0x01       : lg.next[I]
    end

    # Swap current and next grids and return the updated SlowLifeGrid
    lg.grid, lg.next = lg.next, lg.grid
    return lg
end
