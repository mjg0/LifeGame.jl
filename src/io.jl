const MAGIC_BYTES = 0x4546494C



clustertype(::LifeGrid{R,C,H,Tall,Wide}) where {R,C,H,Tall,Wide} = C



# Read and write LifeRules to/from streams
function Base.write(io::IO, ::Type{LifeRule{Rule{B},Rule{S}}}) where {B,S}
    written = write(io, B)
    written += write(io, S)
    return written
end

function LifeRule(io::IO)
    B = read(io, UInt8)
    S = read(io, UInt8)
    return LifeRule{B,S}
end



# Helpers for read and write functions



"""
    write(io::IO, lg::LifeGrid)

Write `lg` to `io` as binary data.

The format is as follows (all in binary, little-endian):

1. 4 magic bytes indicating a non-corrupted `LifeGrid` file/stream
1. The size, 2 64-bit signed integers (height, then width)
1. The rule, 2 bytes representing the birth and survival rules respectively
1. The data grid underpinning `lg`'s cells, a series of 64-bit unsigned integers

All `LifeGrid`s are converted to a cluster size of 64 bits before writing. This won't be a
problem for normal use cases, but could cause a significant allocation if for some reason
you have a large grid with 8-, 16-, or 32-bit clusters.
"""
function Base.write(io::IO, lg::LifeGrid)
    # I/O format is always 64-bit cluster size
    if clustertype(lg) != UInt64
        return write(io, LifeGrid(lg; rule = sprint(show, rule(lg)), CType = UInt64))
    end

    # Write magic bytes
    written = write(io, htol(MAGIC_BYTES))

    # Write size
    written += write(io, htol.(size(lg))...)

    # Write rule
    written += write(io, rule(lg))

    # Write grid, row major, with halos masked to zero for uniformity
    for i in 1:lg.height
        for J in axes(lg.grid, 2)
            masked = lg.grid[i, J] & ~(highbit(UInt64) | lowbit(UInt64))
            written += write(io, htol(masked))
        end
    end

    return written
end



"""
    LifeGrid(io::IO; rule="B3/S23")

Read a `LifeGrid` in from `io` and return it.

See [`write`](@ref) for the binary format of a `LifeGrid`.
"""
function LifeGrid(io::IO; kw...)
    # Ensure the right "magic bytes"
    if ltoh(read(io, typeof(MAGIC_BYTES))) != ltoh(MAGIC_BYTES)
        throw(ArgumentError("Invalid life game grid file: wrong magic bytes"))
    end

    # Read in size
    m = ltoh(read(io, Int64))
    n = ltoh(read(io, Int64))

    # Read in rule
    R = sprint(show, LifeRule(io))

    # Construct the grid to be returned
    lg = LifeGrid(m, n; rule = R, CType = UInt64)

    # Read in cells, row-major
    for i in 1:m
        for J in axes(lg.grid, 2)
            lg.grid[i,J] = ltoh(read(io, UInt64))
        end
    end

    # Return an appropriately sized LifeGrid
    return if size(lg, 2) < nbits(UInt32) - 2
        LifeGrid(lg; rule = R)
    else
        lg
    end
end
