const MAGIC_BYTES = 0x3246494C



clustertype(::LifeGrid{R,C,H,Tall,Wide}) where {R,C,H,Tall,Wide} = C



# Read and write LifeRules to/from streams
function Base.write(io::IO, ::LifeRule{B,S}) where {B,S}
    written = write(io, htol(B))
    written += write(io, htol(S))
    return written
end

function LifeRule(io::IO)
    B = ltoh(read(io, UInt16))
    S = ltoh(read(io, UInt16))
    return LifeRule{B,S}()
end



"""
    write(io::IO, lg::LifeGrid)

Write `lg` to `io` as binary data.

The format is as follows (all in binary, little-endian where applicable):

1. 4 magic bytes indicating a non-corrupted `LifeGrid` file/stream
1. The size, 2 64-bit signed integers (height, then width)
1. The rule, 2 little-endian `UInt16`s representing the birth and survival rules
1. The data grid underpinning `lg`'s cells, packed bits in row-major order

The grid bits are packed with no padding between rows, so a row may begin mid-byte. The grid
is written as 8-bit unsigned integers and is thus endianness-agnostic.
"""
function Base.write(io::IO, lg::LifeGrid)
    # Write magic bytes
    written = write(io, htol(MAGIC_BYTES))

    # Write size
    written += write(io, htol.(size(lg))...)

    # Write rule
    written += write(io, rule(lg))

    # Write grid
    cellspercluster = nbits(eltype(lg.grid)) - 2
    bitwriter(io; bufsize=min(cld(length(lg), 8), 1024^2)) do bw
        @inbounds for i in axes(lg, 1)
            # Write all but the last column
            for j in 1:size(lg.grid, 2)-1
                writebits!(bw, lg.grid[i, j] << 1, cellspercluster)
            end
            # Last column might have less cells
            ncells = size(lg, 2) - (size(lg.grid, 2) - 1) * cellspercluster
            writebits!(bw, lg.grid[i, end] << 1, ncells)
        end
    end
    written += cld(length(lg), 8)

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
        throw(ArgumentError("Invalid life game grid stream: wrong magic bytes"))
    end

    # Read in size
    m = ltoh(read(io, Int64))
    n = ltoh(read(io, Int64))

    # Read in rule
    R = LifeRule(io)

    # Construct the grid to be returned
    lg = LifeGrid(m, n; rule=R, kw...)

    # Read in the grid
    cellspercluster = nbits(eltype(lg.grid)) - 2
    bitreader(io; bufsize=min(cld(length(lg), 8), 1024^2), maxbits=length(lg)) do br
        @inbounds for i in axes(lg, 1)
            # Read all but the last column
            for j in 1:size(lg.grid, 2)-1
                cell = readbits!(br, eltype(lg.grid), cellspercluster)
                lg.grid[i, j] = cell >> 1
                syncclusterhalos!(lg, i, j)
            end
            # Last column might have less cells
            ncells = size(lg, 2) - (size(lg.grid, 2) - 1) * cellspercluster
            cell = readbits!(br, eltype(lg.grid), ncells)
            lg.grid[i, end] = cell >> 1
            syncclusterhalos!(lg, i, lastindex(lg.grid, 2))
        end
    end

    return lg
end
