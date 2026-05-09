export LifePattern



abstract type Mutability end
struct Immutable <: Mutability end
struct Mutable <: Mutability end



struct LifePattern{M<:Mutability} <: AbstractMatrix{Bool}
    height::Int64
    width::Int64
    data::Vector{UInt64} # contains one more than necessary to speed up indexlifepattern

    function LifePattern(m::Integer, n::Integer; mutable = true)
        data = zeros(UInt64, cld(m*n, nbits(UInt64))+1)
        return new{mutable ? Mutable : Immutable}(m, n, data)
    end

    function LifePattern(pattern::AbstractMatrix{<:Number}; mutable = true)
        lp = LifePattern(size(pattern)...)
        lp .= pattern .!= zero(eltype(pattern))
        if mutable
            return lp
        else
            return new{Immutable}(lp.height, lp.width, lp.data)
        end
    end
end



"""
    indexlifepattern(lp::LifePattern, i, j)

Return an index and mask for retrieving a particular bit from `lp.data`.

Returns a ``

Use it thus:

```julia
I, mask = indexlifepattern

# Set to true:
lp.data[I] |=  mask

# Set to false
lp.data[I] &= ~mask
```
"""
Base.@propagate_inbounds @inline function indexlifepattern(lp::LifePattern, i, j)
    T = eltype(lp.data)

    # Flattened index
    ij = (i-1)*lp.width+j

    # lp.data index and appropriate bitshift
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

Base.@propagate_inbounds function Base.setindex!(
    lp::LifePattern{Mutable},
    val::Number,
    i::Integer,
    j::Integer,
)
    I, mask = indexlifepattern(lp, i, j)

    lp.data[I] = ifelse(val==zero(val), lp.data[I] & ~mask, lp.data[I] | mask)

    return lp
end



"""
    bitchunk(lp::LifePattern, i, jrange, ::Type{C})

Return `lp[i:jrange]` as bits packed into the high side of a `C`.
"""
Base.@propagate_inbounds @inline function bitchunk(
    lp::LifePattern,
    i,
    jrange,
    ::Type{C},
) where {C<:Unsigned}
    @boundscheck begin
        checkbounds(lp, i, first(jrange))
        checkbounds(lp, i, last(jrange))
        if last(jrange)-first(jrange)+1 > nbits(C)
            throw(ArgumentError("$jrange is too wide to fit into a $C"))
        end
    end

    data = lp.data
    P = eltype(data)

    # How many bits along are we?
    bit1 = (i-1)*lp.width+first(jrange)
    bit2 = (i-1)*lp.width+last(jrange)

    # How many elements along are we?
    I1 = (bit1-1)÷nbits(P)+1
    I2 = (bit2-1)÷nbits(P)+1
    shift = (bit1-1)%nbits(P)

    # Get all the bits from data[I1:I2] into one P
    x = ifelse(I1==I2, data[I1]<<shift, (data[I1]<<shift) | (data[I2]>>(nbits(P)-shift)))

    # Mask to take off trailing bits
    k = min(bit2-bit1+1, nbits(C))
    mask = ifelse(k==0, zero(P), typemax(P)<<(nbits(P)-k))

    # Return a C, shifted over the appropriate amount
    return C((x & mask) >> (nbits(P)-nbits(C)))
end



"""
    insert(lg::LifeGrid, i::Integer, j::Integer, pattern::AbstractMatrix)
    insert(lg::LifeGrid, I::CartesianIndex{2}, pattern::AbstractMatrix)

Insert `pattern` into `lg`, with the upper left corner at `i, j` or `I`.

Live cells in `pattern` are superimposed over `lg`; dead cells do not overwrite living cells
in `lg`.

If `pattern` is a `LifePattern`, insertion is much more efficient.
"""
Base.@propagate_inbounds function Base.insert!(
    lg::LifeGrid,
    i::Integer,
    j::Integer,
    pattern::AbstractMatrix{<:Number},
)
    for pI in CartesianIndices(pattern)
        I = pI+CartesianIndex(i, j)-oneunit(pI)
        lg[I] = ifelse(pattern[pI]==zero(eltype(pattern)), lg[I], true)
    end

    return lg
end

Base.@propagate_inbounds function Base.insert!(
    lg::LifeGrid{R,C,H},
    i::Integer,
    j::Integer,
    lp::LifePattern,
) where {R,C,H}
    # TODO: @checkbounds

    # Bounds within lg
    I1, J1, lshift = indexlifegrid(lg, i, j)
    I2, J2, rshift = indexlifegrid(lg, ((i, j) .+ size(lp) .- 1)...)

    # Keep track of which column in lp we're on
    pbit = 1

    # Iterate over cells in lg
    for J = J1:J2
        # Where we are, in both lg and lp
        loffset = ifelse(J==J1, lshift, 1)
        roffset = ifelse(J==J2, rshift, nbits(C)-2)
        pJ = pbit:ifelse(J==J2, lastindex(lp, 2), pbit+roffset-loffset)

        # Iterate this row
        for (I, pI) in zip(I1:I2, axes(lp, 1))
            # Update the grid element
            overlay = bitchunk(lp, pI, pJ, C) >> loffset
            lg.grid[I, J] |= overlay

            # Update halos if needed
            hi = i+pI-1
            hj = j+first(pJ)-1
            hidx, hmask = indexhalos(lg, hi, hj)
            if (highbit(C)>>1)&overlay != zero(C)
                lg.halos.currentleft[hidx] |= hmask
            end
            if (lowbit(C)<<1)&overlay != zero(C)
                lg.halos.currentright[hidx] |= hmask
            end
        end

        # Update lp bit start column
        pbit = last(pJ)+1
    end

    return lg
end

Base.@propagate_inbounds function Base.insert!(lg::LifeGrid, I::CartesianIndex{2}, pattern)
    return insert!(lg, Tuple(I)..., pattern)
end
