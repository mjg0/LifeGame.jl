# Types to store a LifeGrid's halo matrices and grid chunk buffers
"""
    Halos{H}

A type to store the halos between the columns of a `LifeGrid`'s `grid`.

There are 4 matrices, each of the same size as `grid`: right and left halos, for both the
current and next iterations.
"""
mutable struct Halos{H}
    currentleft::Matrix{H}
    nextleft::Matrix{H}
    currentright::Matrix{H}
    nextright::Matrix{H}

    Halos(H, m, n) = new{H}(ntuple(_ -> zeros(H, m, n), 4)...)
end



# Bits in `SparseActiveFlags` track which cells in a chunk changed on the last step
const SPARSE_CENTER = UInt16(0x0001)
const SPARSE_TOP = UInt16(0x0002)
const SPARSE_BOTTOM = UInt16(0x0004)
const SPARSE_LEFT = UInt16(0x0008)
const SPARSE_RIGHT = UInt16(0x0010)
const SPARSE_TOPLEFT = UInt16(0x0020)
const SPARSE_TOPRIGHT = UInt16(0x0040)
const SPARSE_BOTTOMLEFT = UInt16(0x0080)
const SPARSE_BOTTOMRIGHT = UInt16(0x0100)



"""
    SparseActiveFlags

Sparse update flags for each chunk of a `LifeGrid`.

`current` is read to decide which chunks need updating; `next` is written with the changes
found during the current step. They are swapped after a sparse step, like [`Halos`](@ref).
`allcurrent` avoids repeatedly filling `current` after dense steps.
"""
mutable struct SparseActiveFlags
    current::Matrix{UInt16}
    next::Matrix{UInt16}
    allcurrent::Bool

    function SparseActiveFlags(m, n)
        current = zeros(UInt16, m + 2, n + 2)
        current[2:m+1, 2:n+1] .= SPARSE_CENTER
        return new(current, zeros(UInt16, m + 2, n + 2), true)
    end
end



"""
    smallestuint(N)

Return the smallest `Unsigned` type that has at least `min(N, 64)` bits.
"""
function smallestuint(N)
    return if N ≤ 8
        UInt8
    elseif N ≤ 16
        UInt16
    elseif N ≤ 32
        UInt32
    else
        UInt64
    end
end
