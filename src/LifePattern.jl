export LifePattern



"""
    LifePattern <: AbstractMatrix{Bool}

A pattern of living and dead cells meant for insertion into a [`LifeGrid`](@ref).

Once a `LifePattern` has been created, it can be inserted into a `LifeGrid` with
[`insert!`](@ref).

Commonly used patterns are available in the [`LifePatterns`](@ref) module.

---

    LifePattern(pattern::Matrix{Bool})
    LifePattern(pattern::BitMatrix)
    LifePattern(pattern::Matrix{<:Number})

Construct a `LifePattern` from a matrix.

Living cells are represented by `true`s, on bits, or non-zero numbers, and dead cells by
`false`s, off bits, or zeros.

The convention is to use a minimal bounding box when creating `LifePattern`s. Oscillators
should have room to transition successfully between states; see for example
[`LifePatterns.pulsar`](@ref).

# Examples

```jldoctest
julia> grid = LifeGrid(5, 8);

julia> pattern = LifePattern([0 1 0
                              1 0 0
                              1 0 1])
3×3 LifePattern:
 0  1  0
 1  0  0
 1  0  1

julia> insert!(grid, 2, 3, pattern)
5×8 LifeGrid:
 0  0  0  0  0  0  0  0
 0  0  0  1  0  0  0  0
 0  0  1  0  0  0  0  0
 0  0  1  0  1  0  0  0
 0  0  0  0  0  0  0  0
```
"""
struct LifePattern <: AbstractMatrix{Bool}
    data::Matrix{UInt8}
    width::Int

    function LifePattern(pattern::Matrix{Bool})
        # Create an array of 8-bit unsigned ints to hold the pattern in its bits
        height, width = size(pattern)
        packedwidth = cld(width, 8)
        data = zeros(UInt8, height, packedwidth)

        # Fill in the grid with bitshifts
        for ij in CartesianIndices(pattern)
            if pattern[ij]
                i, j = Tuple(ij)
                data[i,(j-1)÷8+1] |= 0x80 >> ((j-1)%8)
            end
        end

        # Only the array and width need to be stored, height is part of the array itself
        return new(data, width)
    end

    LifePattern(pattern::BitMatrix) = LifePattern(Matrix{Bool}(pattern))

    LifePattern(pattern::Matrix{T}) where T <: Number = LifePattern(pattern.!=zero(T))
end



# Implement AbstractArray interface for LifePattern
Base.size(pattern::LifePattern) = size(pattern.data, 1), pattern.width

# Get the index of the underlying array given the conceptual index
function indexpattern(i, j)
    I = i
    J = (j-1)÷8+1 # add one for 1-index arrays
    shift = (j-1)%8 # how many bits to shift for the bit in question
    return I, J, shift
end

Base.@propagate_inbounds function Base.getindex(pattern::LifePattern,
                                                i::Integer, j::Integer)
    I, J, shift = indexpattern(i, j)
    return ((pattern.data[I,J] << shift) & 0x80) == 0x80
end

Base.@propagate_inbounds function Base.setindex!(pattern::LifePattern, val::Bool,
                                                 i::Integer, j::Integer)
    I, J, shift = indexpattern(i, j)
    cluster = pattern.data[I,J]
    mask = 0x80 >> shift
    pattern.data[I,J] = ifelse(val,
                               cluster |  mask,
                               cluster & ~mask)
    return val
end

Base.@propagate_inbounds function Base.setindex!(pattern::LifePattern, val::Number,
                                                 i::Integer, j::Integer)
    return pattern[i,j] = val != 0
end



"""
    insert!(lg::LifeGrid, i::Integer, j::Integer, pattern::LifePattern)
    insert!(lg::LifeGrid, I::CartesianIndex{2},   pattern::LifePattern)
    insert!(lg::LifeGrid, i::Integer, j::Integer, pattern::AbstractMatrix)
    insert!(lg::LifeGrid, I::CartesianIndex{2},   pattern::AbstractMatrix)

Insert a pattern into a [`LifeGrid`](@ref) and return the grid.

The upper left corner of the pattern will be inserted at `(i, j)`, with the rest of it to
the right and below.

Using a [`LifePattern`](@ref) is more efficient than arrays of numbers, bools, or bits.

# Examples
```jldoctest
julia> insert!(LifeGrid(3, 5), 2, 2, LifePattern([0 1 0; 1 0 1]))
3×5 LifeGrid:
 0  0  0  0  0
 0  0  1  0  0
 0  1  0  1  0
```
"""
Base.@propagate_inbounds function Base.insert!(lg::LifeGrid, i::Integer, j::Integer,
                                               pattern::LifePattern)
    # We're only interested in the active clusters of lg
    grid = @view lg.grid[begin+1:end-1,begin+1:end-1] # check bounds if the caller wants to

    # Make a bounds check happen
    boundscheckvar = @view lg[i:i+size(pattern, 1)-1,j:j+size(pattern, 2)-1]

    # We need the size of both lg's and pattern's bits
    gridtypebits    = sizeof(eltype(grid        )) * 8
    patterntypebits = sizeof(eltype(pattern.data)) * 8

    # Keep track of how many bits have been inserted into the grid
    insertedbits = 0

    # For each column in pattern.data...
    for y in axes(pattern.data, 2)
        # The number of bits to be inserted may be less than 8 if we're at a cluster border
        bitsthiscolumn = min(patterntypebits, size(pattern, 2)-insertedbits)

        # Keep track of where in pattern.data[i,j] to start if we're at a cluster border
        offset = 0

        # Might need to perform two insertions
        while bitsthiscolumn > 0
            # Figure out which cluster we're on
            I, J, shift = indexlifegrid(i, j+insertedbits)
            I, J = I-1, J-1 # compensate for the view starting at [begin+1,begin+1]

            # How many bits need to be inserted into this cluster?
            bitstoboundary = CELLS_PER_CLUSTER - shift + 1
            toinsert = min(bitsthiscolumn, bitstoboundary)

            # Do the same thing for each row
            @simd for x in axes(pattern.data, 1)
                # Which bits of the cluster are up for replacement?
                mask = (typemax(CLUSTER_TYPE) << (gridtypebits-toinsert)) >> shift
                # Get the actual values to be inserted
                insert = ((CLUSTER_TYPE(pattern.data[x,y] << offset)
                           << (gridtypebits-patterntypebits)) >> shift) & mask
                # Insert the appropriate portion of the pattern
                grid[I+x-1, J] = (grid[I+x-1, J] & ~mask) | insert
            end

            # Update iteration parameters
            offset = toinsert
            bitsthiscolumn -= toinsert
            insertedbits += toinsert
        end
    end

    return lg
end



# Implement AbstractArray interface for LifePattern
Base.@propagate_inbounds function Base.insert!(lg::LifeGrid, I::CartesianIndex{2},
                                               pattern::LifePattern)
    return insert!(lg, Tuple(I)..., pattern)
end

Base.@propagate_inbounds function Base.insert!(lg::LifeGrid, i::Integer, j::Integer,
                                               pattern::AbstractMatrix)
    m, n = size(pattern)
    lg[i:i+m-1,j:j+n-1] .= pattern
    return lg
end

Base.@propagate_inbounds function Base.insert!(lg::LifeGrid, I::CartesianIndex{2},
                                               pattern::AbstractMatrix)
    return insert!(lg, Tuple(I)..., pattern)
end