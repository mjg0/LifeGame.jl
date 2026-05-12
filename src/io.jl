const MAGIC_BYTES = 0x4546494C



clustertype(::LifeGrid{R,C,H,Tall,Wide}) where {R,C,H,Tall,Wide} = C



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



function Base.write(io::IO, lg::LifeGrid)
    # I/O format is always 64-bit cluster size
    if clustertype(lg) != UInt64
        return write(io, LifeGrid(lg; rule = sprint(show, rule(lg)), CType = UInt64))
    end

    # Write magic bytes
    written = write(io, MAGIC_BYTES)

    # Write size
    written += write(io, htol.(size(lg))...)

    # Write rule
    written += write(io, rule(lg))

    # Write grid, with halos masked to zero for uniformity
    for J in axes(lg.grid, 2)
        for i in 1:lg.height
            masked = lg.grid[i, J] & ~(highbit(UInt64) | lowbit(UInt64))
            written += write(io, htol(masked))
        end
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
    R = sprint(show, LifeRule(io))

    # Construct the grid to be returned
    lg = LifeGrid(m, n; rule = R, CType = UInt64)

    # Read in cells
    for J in axes(lg.grid, 2)
        for i in 1:m
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
