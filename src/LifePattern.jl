export LifePattern



struct LifePattern <: AbstractMatrix{Bool}
    height::Int64
    width::Int64
    data::Vector{UInt64}

    function LifePattern(m::Integer, n::Integer)
        data = zeros(UInt64, cld(m*n, nbits(UInt64))+1) # store a one past the end
        return new(m, n, data)
    end

    function LifePattern(pattern::AbstractMatrix{<:Number})
        lp = LifePattern(size(pattern)...)
        lp .= pattern.!=zero(eltype(pattern))
        return lp
    end
end



# Return the index and mask for pattern.data given conceptual index (i, j)
Base.@propagate_inbounds @inline function indexlifepattern(lp::LifePattern, i, j, T::Type{<:Unsigned}=UInt64)
    ij = (i-1)*lp.width+j

    I = (ij-1)÷nbits(T)+1
    k = (ij-1)%nbits(T)

    mask = highbit(T) >> k

    return I, mask
end



# AbstractArray interface for LifePattern
Base.size(lp::LifePattern) = lp.height, lp.width

Base.@propagate_inbounds function Base.getindex(lp::LifePattern, i::Integer, j::Integer)
    I, mask = indexlifepattern(lp, i, j)

    return mask & lp.data[I] != zero(UInt64)
end

Base.@propagate_inbounds function Base.setindex!(lp::LifePattern, val::Number, i::Integer, j::Integer)
    I, mask = indexlifepattern(lp, i, j)

    lp.data[I] = ifelse(val==zero(val),
                        lp.data[I],
                        lp.data[I] | mask)

    return lp
end



Base.@propagate_inbounds @inline function bitchunk(lp::LifePattern, i, jrange, ::Type{C}) where C
    # TODO: @checkinbounds if jrange is small enough to fit in C
    # TODO: assert nbits(C)≤nbits(UInt64)

    # Which bits in the lp.data vector are we working with?
    bitrange = (i-1)*lp.width.+jrange

    # Convenience array to see the data as type C
    data = reinterpret(C, lp.data)

    # Which indices in data we're working with
    I1 = first(bitrange)÷nbits(C)+1
    I2 = last( bitrange)÷nbits(C)+1

    # How many bits to drop from the front of data[I1]
    shift = (first(bitrange)-1)÷nbits(C)+1

    # Get all the bits from data[I1:I2] into one C
    x = ifelse(I1==I2, # This necessitates the extra cell in lp.data
               data[I1]<<shift,
               data[I1]<<shift | data[I2]>>(nbits(C)-shift))

    # Mask off the trailing bits
    k = clamp(last(bitrange)-first(bitrange)+1, 0, nbits(C))
    mask = ifelse(k==0,
                  zero(C),
                  typemax(C)<<(nbits(C)-k))

    return x & mask
end




# Speed to beat: about 20x faster than naive insert with a sparse matrix, 1000x1000
function Base.insert!(lg::LifeGrid, i::Integer, j::Integer, pattern::AbstractMatrix{<:Number})
    for pI in CartesianIndices(pattern)
        I = pI+CartesianIndex(i, j)-one(pI)
        lg[I] = ifelse(pattern[pI]==zero(eltype(pI)),
                       lg[I],
                       true)
    end
end

function Base.insert!(lg::LifeGrid{R, C, H}, i::Integer, j::Integer, lp::LifePattern) where {R, C, H}
    # TODO: @checkbounds

    # Bounds within lg
    I1, J1, lshift = indexlifegrid(lg, i, j)
    I2, J2, rshift = indexlifegrid(lg, i+size(lp, 1)-1, j+size(lp, 2)-1)

    # Keep track of which column in lp we're on
    pbit = 1

    # Iterate over cells in lg
    for J in J1:J2
        # Where we are, in both lg and lp
        loffset =  ifelse(J==J1, lshift,                1)
        roffset =  ifelse(J==J2, rshift,                nbits(C)-2)
        pJ = pbit:ifelse(J==J2, lastindex(lp, 2), pbit+roffset-loffset+1)

        # Iterate this row
        for (I, pI) in zip(I1:I2, axes(lp, 1))
            # Overlay the appropriate bits
            overlay = bitchunk(lp, pI, pJ, C) >> loffset
            lg.grid[I,J] |= overlay

            # Update halos
            hidx, hmask = indexhalos(lg, I, J)
            if (highbit(C) >> 1) & overlay != zero(C)
                lg.lefthalos[1][hidx] |= hmask
            end
            if (lowbit( C) << 1) & overlay != zero(C)
                lg.righthalos[1][hidx] |= hmask
            end
        end

        # Update lp bit start column
        pbit = last(pJ)+1
    end

    return lg
end