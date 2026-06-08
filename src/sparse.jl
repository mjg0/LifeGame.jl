"""
    sparsechangedbits(delta, i, H)

Return sparse-update flags for the `i`th cluster of a chunk.

`delta` has on bits where cells changed in this cluster. Any change sets the center flag;
changes in the first or last real cell set horizontal flags; changes in the first or last
cluster set vertical flags. Neighboring chunks read these flags on the next sparse step.
"""
Base.@propagate_inbounds @inline function sparsechangedbits(delta::C) where {C}
    if delta == zero(C)
        return zero(UInt16)
    end

    leftchanged = delta & (highbit(C) >> 1) != zero(C)
    rightchanged = delta & (lowbit(C) << 1) != zero(C)
    return SPARSE_CENTER |
        ifelse(leftchanged, SPARSE_LEFT, zero(UInt16)) |
        ifelse(rightchanged, SPARSE_RIGHT, zero(UInt16))
end

Base.@propagate_inbounds @inline function sparsechangedbits(delta::C, i, ::Type{H}) where {C,H}
    bits = sparsechangedbits(delta)
    bits == zero(UInt16) && return bits

    if i == 1
        bits |= SPARSE_TOP
        bits |= ifelse(bits & SPARSE_LEFT != zero(UInt16), SPARSE_TOPLEFT, zero(UInt16))
        bits |= ifelse(bits & SPARSE_RIGHT != zero(UInt16), SPARSE_TOPRIGHT, zero(UInt16))
    elseif i == nbits(H)
        bits |= SPARSE_BOTTOM
        bits |= ifelse(bits & SPARSE_LEFT != zero(UInt16), SPARSE_BOTTOMLEFT, zero(UInt16))
        bits |= ifelse(bits & SPARSE_RIGHT != zero(UInt16), SPARSE_BOTTOMRIGHT, zero(UInt16))
    end

    return bits
end

Base.@propagate_inbounds @inline function sparsechangedbits(
    delta::C,
    topdelta::C,
    bottomdelta::C,
) where {C}
    delta == zero(C) && return zero(UInt16)

    leftbit = highbit(C) >> 1
    rightbit = lowbit(C) << 1

    bits = SPARSE_CENTER
    bits |= ifelse(delta & leftbit != zero(C), SPARSE_LEFT, zero(UInt16))
    bits |= ifelse(delta & rightbit != zero(C), SPARSE_RIGHT, zero(UInt16))
    bits |= ifelse(topdelta != zero(C), SPARSE_TOP, zero(UInt16))
    bits |= ifelse(topdelta & leftbit != zero(C), SPARSE_TOPLEFT, zero(UInt16))
    bits |= ifelse(topdelta & rightbit != zero(C), SPARSE_TOPRIGHT, zero(UInt16))
    bits |= ifelse(bottomdelta != zero(C), SPARSE_BOTTOM, zero(UInt16))
    bits |= ifelse(
        bottomdelta & leftbit != zero(C),
        SPARSE_BOTTOMLEFT,
        zero(UInt16),
    )
    bits |= ifelse(
        bottomdelta & rightbit != zero(C),
        SPARSE_BOTTOMRIGHT,
        zero(UInt16),
    )

    return bits
end



"""
    chunkactive(changed, I, J)

Return whether chunk `(I, J)` needs updating from the sparse flags in `changed`.
"""
Base.@propagate_inbounds @inline function chunkactive(changed::Matrix{UInt16}, I, J)
    I += 1
    J += 1
    changed[I, J] & SPARSE_CENTER != zero(UInt16) && return true

    return (
        changed[I - 1, J] & SPARSE_BOTTOM |
        changed[I + 1, J] & SPARSE_TOP |
        changed[I, J - 1] & SPARSE_RIGHT |
        changed[I, J + 1] & SPARSE_LEFT |
        changed[I - 1, J - 1] & SPARSE_BOTTOMRIGHT |
        changed[I - 1, J + 1] & SPARSE_BOTTOMLEFT |
        changed[I + 1, J - 1] & SPARSE_TOPRIGHT |
        changed[I + 1, J + 1] & SPARSE_TOPLEFT
    ) != zero(UInt16)
end



"""
    marksparsechanged!(lg::LifeGrid, I, J, delta)

Update the sparse activity tracking bits for the `(I, J)`th chunk of `lg` given `delta`.
"""
Base.@propagate_inbounds @inline function marksparsechanged!(
    lg::LifeGrid{R,C,H},
    I::Integer,
    J::Integer,
    delta::C,
) where {R,C,H}
    if lg.changed.allcurrent || delta == zero(C)
        return nothing
    end

    hI = (I - 1) ÷ nbits(H) + 1
    i = (I - 1) % nbits(H) + 1
    lg.changed.current[hI + 1, J + 1] |= sparsechangedbits(delta, i, H)

    return nothing
end
