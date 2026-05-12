


const MAGIC_BYTES = 0x4C494645

clustertype(::LifeGrid{R,C,H,Tall,Wide}) where {R,C,H,Tall,Wide} = C



function Base.write(io::IO, ::LifeRule{B,S}) where {B,S}
    write(io, B)
    write(io, S)
end

function LifeRule(io::IO)
    B = read(io, UInt8)
    S = read(io, UInt8)
    return LifeRule{B,S}()
end



function Base.write(io::IO, lg::LifeGrid)
    # I/O format is always 64-bit cluster size
    if clustertype(lg) != UInt64
        return write(io, LifeGrid(lg; rule=rule(lg), CType=UInt64))
    end

    # Write magic bytes
    written = write(io, MAGIC_BYTES)

    # Write size
    written += write(io, htol.(size(lg))...)

    # Write rule
    written += write(io, rule(lg))

    # Write grid, with halos masked to zero for uniformity
    for cluster in lg.grid
        masked = cluster & ~(highbit(UInt64) | lowbit(UInt64))
        written += write(io, htol(masked))
    end

    return written
end



function LifeGrid(io::IO; kw...)
    # Ensure the right "magic bytes"
    if read(io, typeof(MAGIC_BYTES)) != MAGIC_BYTES
        throw(ArgumentError("Invalid life game grid file: wrong magic bytes"))
    end

    # Read in size
    m = ltoh(read(io, Int64))
    n = ltoh(read(io, Int64))

    # Read in rule
    R = rule(read(io, String))

    # Construct the grid to be returned
    lg = LifeGrid(m, n; rule=R, CType=UInt64)

    # Read in cells
    read!(io, lg.grid)
    map!(ltoh, grid, grid)

    # Return an appropriately sized LifeGrid
    return if size(lg, 2) < nbits(UInt32) - 2
        LifeGrid(lg; rule=R)
    else
        lg
    end
end