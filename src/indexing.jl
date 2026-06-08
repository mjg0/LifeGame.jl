"""
    indexlifegrid(lg::LifeGrid, i, j)

Return the coordinates and a shift translating index `(i, j)` to a cell in `lg.grid`.

Coordinates are two indices and the shift is a number up to `8*sizeof(eltype(lg.grid))`.

For coordinates and shift obtained thus:

```julia
I, J, shift = indexlifegrid(lg, i, j)
```

...`lg[i,j]` is true only if `(lg.grid[I,J] << shift` has its highest bit on.
"""
Base.@propagate_inbounds @inline function indexlifegrid(
    ::LifeGrid{R,C,H},
    i,
    j,
) where {R,C,H}
    I = i
    J = (j - 1) ÷ (nbits(C) - 2) + 1
    shift = (j - 1) % (nbits(C) - 2) + 1

    return I, J, shift
end



"""
    gridchunk(lg::LifeGrid, I, J)

Return a view of the `(I, J)`th chunk of `lg`.

A "chunk" is a section of a column as many elements long as `lg`'s halo type has bits. This
is usually 64 clusters for large grids, but could be as few as 8.
"""
Base.@propagate_inbounds @inline function gridchunk(
    lg::LifeGrid{R,C,H,Tall,Wide},
    I,
    J,
) where {R,C,H,Tall,Wide}
    i = (I - 1) * nbits(H) + 1
    return view(lg.grid, i:(i+nbits(H)-1), J)
end



"""
    indexhalos(lg::LifeGrid, i, j)

Helper for `getindex(::LifeGrid, args...)` and `setindex!(::LifeGrid, args...)`.

Return an index and mask for indexing the halos in `lg`, given conceptual index `(i, j)`.
"""
Base.@propagate_inbounds @inline function indexhalos(
    lg::LifeGrid{R,C,H},
    i::Integer,
    j::Integer,
) where {R,C,H}
    I, J, _ = indexlifegrid(lg, i, j)

    # The index in the halo array
    hidx = CartesianIndex((I - 1) ÷ nbits(H) + 1, J)

    # A single-bit mask with the appropriate shift
    k = (I - 1) % nbits(H) + 1
    mask = one(H) << (nbits(H) - k)

    return hidx, mask
end



"""
    syncclusterhalos!(lg::LifeGrid, I, J)

Update `lg`'s current halo matrices for the `(I, J)`th cluster of `lg.grid`.

Helper for `setindex!(::LifeGrid, args...)`.

This should be called after directly mutating `lg.grid[I, J]` without using `setindex!`. It
copies the first and last real cell bits of the cluster into `lg.halos.currentleft` and
`lg.halos.currentright`.
"""
Base.@propagate_inbounds function syncclusterhalos!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    I::Integer,
    J::Integer,
) where {R,C,H,Tall,Wide}
    hidx = CartesianIndex((I - 1) ÷ nbits(H) + 1, J)
    k = (I - 1) % nbits(H) + 1
    hmask = one(H) << (nbits(H) - k)

    cluster = lg.grid[I, J]

    leftbit = (cluster >> (nbits(C) - 2)) & one(C)
    rightbit = (cluster >> 1) & one(C)

    lg.halos.currentleft[hidx] = ifelse(
        leftbit == one(C),
        lg.halos.currentleft[hidx] | hmask,
        lg.halos.currentleft[hidx] & ~hmask,
    )

    lg.halos.currentright[hidx] = ifelse(
        rightbit == one(C),
        lg.halos.currentright[hidx] | hmask,
        lg.halos.currentright[hidx] & ~hmask,
    )

    return nothing
end



# AbstractArray interface for LifeGrid
Base.size(lg::LifeGrid) = lg.height, lg.width

Base.@propagate_inbounds function Base.getindex(
    lg::LifeGrid{R,C,H},
    i::Integer,
    j::Integer,
) where {R,C,H}
    I, J, shift = indexlifegrid(lg, i, j)

    return ((lg.grid[I, J] << shift) & highbit(C)) == highbit(C)
end

Base.@propagate_inbounds function Base.setindex!(
    lg::LifeGrid{R,C,H},
    val::Number,
    i::Integer,
    j::Integer,
) where {R,C,H}
    I, J, shift = indexlifegrid(lg, i, j)

    # Update the cluster
    cellmask = highbit(C) >> shift
    cluster = lg.grid[I, J]
    newcluster = ifelse(val != zero(val), cluster | cellmask, cluster & ~cellmask)
    lg.grid[I, J] = newcluster
    marksparsechanged!(lg, I, J, cluster ⊻ newcluster)

    # Update halos if needed
    if shift == 1 || shift == nbits(C) - 2
        syncclusterhalos!(lg, I, J)
    end

    return val
end
