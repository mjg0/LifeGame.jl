"""
    insert(lg::LifeGrid, i::Integer, j::Integer, pattern::AbstractMatrix)
    insert(lg::LifeGrid, I::CartesianIndex{2}, pattern::AbstractMatrix)

Insert `pattern` into `lg`, with its upper left corner at `lg[i, j]` or `lg[I]`.

Only living cells in `pattern` overwrite cells in `lg`.

Insertion is much more efficient if `pattern` is a `LifePattern`.
"""
Base.@propagate_inbounds function Base.insert!(
    lg::LifeGrid,
    i::Integer,
    j::Integer,
    pattern::AbstractMatrix{<:Number},
)
    for pI in CartesianIndices(pattern)
        I = pI + CartesianIndex(i, j) - oneunit(pI)
        lg[I] = ifelse(pattern[pI] == zero(eltype(pattern)), lg[I], true)
    end

    return lg
end

Base.@propagate_inbounds function Base.insert!(
    lg::LifeGrid{R,C,H},
    i::Integer,
    j::Integer,
    lp::LifePattern,
) where {R,C,H}
    @boundscheck checkbounds(lg, ((i, j) .+ size(lp) .- 1)...)

    # Bounds within lg
    I1, J1, lshift = indexlifegrid(lg, i, j)
    I2, J2, rshift = indexlifegrid(lg, ((i, j) .+ size(lp) .- 1)...)

    # Keep track of which column in lp we're on
    pbit = 1

    # Iterate over cells in lg
    for J = J1:J2
        # Where we are, in both lg and lp
        loffset = ifelse(J == J1, lshift, 1)
        roffset = ifelse(J == J2, rshift, nbits(C) - 2)
        pJ = pbit:ifelse(J==J2, lastindex(lp, 2), pbit+roffset-loffset)

        # Iterate this row
        for (I, pI) in zip(I1:I2, axes(lp, 1))
            # Update the grid element
            overlay = bitchunk(lp, pI, pJ, C) >> loffset
            delta = overlay & ~lg.grid[I, J]
            lg.grid[I, J] |= overlay
            marksparsechanged!(lg, I, J, delta)

            # Update halos if needed
            hi = i + pI - 1
            hj = j + first(pJ) - 1
            hidx, hmask = indexhalos(lg, hi, hj)
            if (highbit(C) >> 1) & overlay != zero(C)
                lg.halos.currentleft[hidx] |= hmask
            end
            if (lowbit(C) << 1) & overlay != zero(C)
                lg.halos.currentright[hidx] |= hmask
            end
        end

        # Update lp bit start column
        pbit = last(pJ) + 1
    end

    return lg
end

Base.@propagate_inbounds function Base.insert!(lg::LifeGrid, I::CartesianIndex{2}, pattern)
    return insert!(lg, Tuple(I)..., pattern)
end
