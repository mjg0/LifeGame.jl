export LifePattern



struct LifePattern <: AbstractMatrix{Bool}
    height::Int64
    width::Int64
    data::Vector{UInt64} # contains one more than necessary to speed up indexlifepattern

    function LifePattern(m::Integer, n::Integer)
        data = zeros(UInt64, cld(m*n, nbits(UInt64))+1)
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

Base.@propagate_inbounds function Base.setindex!(lp::LifePattern, val::Number, i::Integer, j::Integer)
    I, mask = indexlifepattern(lp, i, j)

    lp.data[I] = ifelse(val==zero(val),
                        lp.data[I] & ~mask,
                        lp.data[I] |  mask)

    return lp
end





#Base.@propagate_inbounds @inline function bitchunk(lp::LifePattern, i, jrange, ::Type{C}) where C
#    # TODO: @checkinbounds if jrange is small enough to fit in C
#    # TODO: assert nbits(C)≤nbits(UInt64)
#
#    # Which bits in the lp.data vector are we working with?
#    bitrange = (i-1)*lp.width.+jrange
#
#    # Convenience array to see the data as type C
#    data = reinterpret(C, lp.data)
#
#    # Which indices in data we're working with
#    I1 = first(bitrange)÷nbits(C)+1
#    I2 = last( bitrange)÷nbits(C)+1
#
#    # How many bits to drop from the front of data[I1]
#    shift = (first(bitrange)-1)%nbits(C)
#
#    # Get all the bits from data[I1:I2] into one C
#    println(data[I1]|>bitstring)
#    x = ifelse(I1==I2, # This necessitates the extra cell in lp.data
#               data[I1]<<shift,
#               data[I1]<<shift | data[I2]>>(nbits(C)-shift))
#
#    # Mask off the trailing bits
#    k = clamp(last(bitrange)-first(bitrange)+1, 0, nbits(C))
#    mask = ifelse(k==0,
#                  zero(C),
#                  typemax(C)<<(nbits(C)-k))
#
#    println(I1:I2)
#    println(bitrange, ' ', shift, ' ', bitstring(x), ' ', bitstring(mask), ' ', bitstring(x&mask))
#    return x & mask
#end




# Speed to beat: about 20x faster than naive insert with a sparse matrix, 1000x1000
function Base.insert!(lg::LifeGrid, i::Integer, j::Integer, pattern::AbstractMatrix{<:Number})
    for pI in CartesianIndices(pattern)
        I = pI+CartesianIndex(i, j)-oneunit(pI)
        lg[I] = ifelse(pattern[pI]==zero(eltype(pI)),
                       lg[I],
                       true)
    end

    return lg
end

# Use this if each LifePattern row is padded to a UInt64 boundary.
@inline bitstride(lp::LifePattern) = nbits(UInt64) * cld(lp.width, nbits(UInt64))

Base.@propagate_inbounds @inline function bitchunk(
    lp::LifePattern,
    i,
    jrange,
    ::Type{C},
) where {C<:Unsigned}
    Cbits = nbits(C)
    Wbits = nbits(UInt64)

    k = last(jrange) - first(jrange) + 1
    k = clamp(k, 0, Cbits)
    k == 0 && return zero(C)

    # Logical bit indices, 1-based, MSB-first within each UInt64.
    #
    # If your LifePattern rows are packed back-to-back with no row padding,
    # replace `bitstride(lp)` with `lp.width`.
    bit1 = (i - 1) * bitstride(lp) + first(jrange)

    w0, r = divrem(bit1 - 1, Wbits)
    w = w0 + 1

    # Pull up to 64 logical MSB-first bits into the top of x64.
    @inbounds x64 = lp.data[w] << r

    if r != 0
        @inbounds x64 |= lp.data[w + 1] >> (Wbits - r)
    end

    # Project the high Cbits into C.
    x = if Cbits == Wbits
        C(x64)
    else
        C(x64 >> (Wbits - Cbits))
    end

    # Keep only the requested logical bits, left-aligned.
    mask = typemax(C) << (Cbits - k)

    return x & mask
end