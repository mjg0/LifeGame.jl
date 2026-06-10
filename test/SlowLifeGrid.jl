using LifeGame # need to overload LifeGame.step!



# Basic life grid implementation against which to test LifeGrid
mutable struct SlowLifeGrid <: AbstractMatrix{Bool}
    grid::Matrix{Bool} # the grid itself
    next::Matrix{Bool} # area to store intermediate results
    birthsums::Vector{Int}    # list of neighbor sums that lead to cell birth
    survivalsums::Vector{Int} # list of neighbor sums that allow cell survival

    function SlowLifeGrid(height, width; rule="B3/S23")
        rulematch = match(r"^B(\d*)/S(\d*)$", rule)

        birth = [parse(Int, c) for c in first(rulematch.captures)]
        survival = [parse(Int, c) for c in last(rulematch.captures)]

        return new(zeros(height, width), zeros(height, width), birth, survival)
    end

    SlowLifeGrid(grid; kw...) = SlowLifeGrid(size(grid)...; kw...) .= grid
end



# Implement AbstractArray interface for SlowLifeGrid
Base.size(lg::SlowLifeGrid) = size(lg.grid)

Base.@propagate_inbounds function Base.getindex(lg::SlowLifeGrid, x...)
    return getindex(lg.grid, x...)
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
        I1 = max(first(region), I - oneunit(I))
        I2 = min(last(region),  I + oneunit(I))
        neighborsum = sum(@view lg.grid[I1:I2]) - lg.grid[I]

        survival = neighborsum in lg.survivalsums
        birth = neighborsum in lg.birthsums

        lg.next[I] = lg.grid[I] ? survival : birth
    end

    # Swap current and next grids and return the updated SlowLifeGrid
    lg.grid, lg.next = lg.next, lg.grid
    return lg
end



# Insert a pattern into a SlowLifeGrid
function Base.insert!(lg::SlowLifeGrid, I::CartesianIndex{2}, pattern)
    for pI in CartesianIndices(pattern)
        idx = pI + I - oneunit(pI)
        lg[idx] = ifelse(pattern[pI] == zero(eltype(pattern)), lg[idx], true)
    end

    return lg
end

function Base.insert!(lg::SlowLifeGrid, i::Integer, j::Integer, pattern)
    return insert!(lg, CartesianIndex((i, j)), pattern)
end
