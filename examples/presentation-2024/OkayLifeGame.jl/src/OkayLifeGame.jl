module OkayLifeGame

export OkayLifeGrid, step!



# Life grid with an array to hold sums of neighboring cells
struct OkayLifeGrid <: AbstractMatrix{Bool}
    grid::Matrix{UInt8}
    sums::Matrix{UInt8}

    OkayLifeGrid(height, width) = new(zeros(Bool, height+2, width+2),
                                      zeros(Bool, height,   width  ))

    OkayLifeGrid(grid) = OkayLifeGrid(size(grid)...) .= grid
end



# Implement AbstractArray interface for OkayLifeGrid
Base.size(lg::OkayLifeGrid) = size(lg.grid).-2

Base.@propagate_inbounds function Base.getindex( lg::OkayLifeGrid, i::Integer, j::Integer)
    return !iszero(lg.grid[i+1,j+1])
end

Base.@propagate_inbounds function Base.setindex!(lg::OkayLifeGrid, val::Bool,
                                                 i::Integer, j::Integer)
    return lg.grid[i+1,j+1] = reinterpret(UInt8, val)
end

Base.@propagate_inbounds function Base.setindex!(lg::OkayLifeGrid, val::Number,
                                                 i::Integer, j::Integer)
    return lg.grid[i+1,j+1] = !iszero(val)
end



# Update an OkayLifeGrid
function step!(lg::OkayLifeGrid)
    # Views of "central" and neighboring cells
    ((upperleft, upper, upperright),
     (     left, center,     right),
     (lowerleft, lower, lowerright)) = ntuple(J->begin
        ntuple(I->begin
            offsets = (I-2, J-2)
            region = ntuple(N->firstindex(lg.grid, N)+1+offsets[N]
                               :lastindex(lg.grid, N)-1+offsets[N], 2)
            return view(lg.grid, region...)
        end, 3)
    end, 3)

    # Calculate neighbor sums and store them in lg.sums
    @inbounds @simd for I in eachindex(center)
        lg.sums[I] = center[I] << 4 # Store this cell's state in the 5th bit
        lg.sums[I] += (upperleft[I] + upper[I] + upperright[I] +
                            left[I] +                 right[I] +
                       lowerleft[I] + lower[I] + lowerright[I])
    end

    # Update center cells based on lg.sums
    @inbounds @simd for I in eachindex(center)
        alive = lg.sums[I] == 0x13 || lg.sums[I] == 0x12 || lg.sums[I] == 0x03
        center[I] = ifelse(alive, 1, 0)
    end

    # Return the grid
    return lg
end



end # module
