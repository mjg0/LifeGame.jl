# Helper functions
@inline realcells(cluster::C) where {C} = cluster & ~(lowbit(C) | highbit(C))

Base.@propagate_inbounds @inline function incominghalos(
    lg::LifeGrid{R,C,H,Tall,Wide},
    I,
    J,
) where {R,C,H,Tall,Wide}
    lhalo = J == 1 ? zero(H) : lg.halos.currentright[I, J-1]
    rhalo = J == size(lg.grid, 2) ? zero(H) : lg.halos.currentleft[I, J+1]
    return lhalo, rhalo
end

Base.@propagate_inbounds @inline function horizontalhalobits(
    ::LifeGrid{R,C,H,Tall,false},
    lhalo::H,
    rhalo::H,
    k,
) where {R,C,H,Tall}
    return zero(C)
end

Base.@propagate_inbounds @inline function horizontalhalobits(
    ::LifeGrid{R,C,H,Tall,true},
    lhalo::H,
    rhalo::H,
    k,
) where {R,C,H,Tall}
    lbit = C(lhalo >> (nbits(H) - k) & one(H)) << (nbits(C) - 1)
    rbit = C(rhalo >> (nbits(H) - k) & one(H))
    return lbit | rbit
end

Base.@propagate_inbounds @inline function addhalos(lg, cluster, lhalo, rhalo, k)
    return realcells(cluster) | horizontalhalobits(lg, lhalo, rhalo, k)
end



"""
    updatebuffers!(buffer::Vector, lg::LifeGrid, I, J)

Update the halo bits of the `(I, J)`th chunk of `lg`, storing the result in `buffer`.

`buffer` should have two elements more than the halo type of `lg` has bits. The first and
last elements should store clusters from the previous and next buffers, respectively.
"""
Base.@propagate_inbounds @inline function updatebuffers!(
    buffer::Vector{C},
    lg::LifeGrid{R,C,H,Tall,Wide},
    I::Integer,
    J::Integer,
) where {R,C,H,Tall,Wide}
    @boundscheck if length(buffer) != nbits(H) + 2
        throw(ArgumentError("buffer should have nbits(H)+2 elements"))
    end

    chunk = gridchunk(lg, I, J)
    lhalo, rhalo = incominghalos(lg, I, J)

    @simd for k = 1:nbits(H)
        buffer[k+1] = addhalos(lg, chunk[k], lhalo, rhalo, k)
    end

    return nothing
end



"""
    updategridchunk!(lg::LifeGrid, buffer::Vector, I, J)

Update the clusters of `buffer` and store the results in the `(I, J)`th chunk of `lg`.

Since `updatedcluster` requires clusters above and below, `buffer` should have two more
elements than the halo type of `lg` has bits. The first and last elements must be
initialized by the caller.
"""
Base.@propagate_inbounds @inline function updategridchunk!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    buffer::Vector{C},
    I::Integer,
    J::Integer,
) where {R,C,H,Tall,Wide}
    return updategridchunk!(lg, buffer, I, J, sparse)
end

Base.@propagate_inbounds @inline function updategridchunk!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    buffer::Vector{C},
    I::Integer,
    J::Integer,
    ::Type{Sp},
) where {R,C,H,Tall,Wide,Sp<:Sparsity}
    chunk = gridchunk(lg, I, J)

    # How many bits to roll off the last column
    trailingbits = ifelse(J == size(lg.grid, 2), mod(-size(lg, 2), nbits(C) - 2) + 1, 0)

    # Sparse chunk activity tracking
    centermask = ~(lowbit(C) | highbit(C))
    delta = zero(C)
    topdelta = zero(C)
    bottomdelta = zero(C)

    # Update each cluster in chunk
    @simd for i = 1:nbits(H)
        # Calculate the updated cluster
        updated = updatedcluster(buffer[i], buffer[i+1], buffer[i+2], R())

        # Keep padding bits out of sparse change detection
        updated = (updated >> trailingbits) << trailingbits

        # Sparse active chunk tracker update
        if Sp === sparse
            rowdelta = (buffer[i+1] ⊻ updated) & centermask
            delta |= rowdelta
            topdelta |= ifelse(i == 1, rowdelta, zero(C))
            bottomdelta |= ifelse(i == nbits(H), rowdelta, zero(C))
        end

        chunk[i] = updated
    end

    # Zero trailing rows if needed
    if I == size(lg.halos.currentleft, 1)
        lg.grid[size(lg, 1)+1:end, J] .= zero(C)
    end

    # Sparse bookkeeping
    if Sp === sparse
        lg.changed.next[I + 1, J + 1] = sparsechangedbits(delta, topdelta, bottomdelta)
    elseif !lg.changed.allcurrent
        lg.changed.current[I + 1, J + 1] = SPARSE_CENTER
    end

    return nothing
end



"""
    updatehalos!(lg::LifeGrid, I, J)

Read the halo bits from the `(I, J)`th chunk of `lg` and put them in `lg.halos.next*`.

The halos are not switched in this function; they should be swapped manually at the end of
a grid update.
"""
Base.@propagate_inbounds @inline function updatehalos!(
    lg::LifeGrid{R,C,H,Tall,true},
    I::Integer,
    J::Integer,
) where {R,C,H,Tall}
    chunk = gridchunk(lg, I, J)

    lhalo = zero(H)
    rhalo = zero(H)

    for k = 1:nbits(H)
        # Get this cluster's halos
        lshift = nbits(C) - 2
        rshift = 1
        lbit = H((chunk[k] >> lshift) & one(C))
        rbit = H((chunk[k] >> rshift) & one(C))

        # Update lhalo and rhalo at this iteration's bit
        lhalo |= lbit << (nbits(H) - k)
        rhalo |= rbit << (nbits(H) - k)
    end

    lg.halos.nextleft[I, J] = lhalo
    lg.halos.nextright[I, J] = rhalo

    return nothing
end

# No need for halo updates on single-column grids
updatehalos!(::LifeGrid{R,C,H,Tall,false}, args...) where {R,C,H,Tall} = nothing



"""
    haloedcluster(lg::LifeGrid, I, J, k)

Return one cluster with current horizontal halo bits applied.

This is used by sparse updates to fill the buffer row above or below an active chunk when
the neighboring chunk was skipped.
"""
Base.@propagate_inbounds @inline function haloedcluster(
    lg::LifeGrid{R,C,H,Tall,Wide},
    I::Integer,
    J::Integer,
    k::Integer,
) where {R,C,H,Tall,Wide}
    if I < 1 || I > size(lg.halos.currentleft, 1)
        return zero(C)
    end

    cluster = lg.grid[(I - 1) * nbits(H) + k, J]
    lhalo, rhalo = incominghalos(lg, I, J)
    return addhalos(lg, cluster, lhalo, rhalo, k)
end



"""
    updatebufferedchunk!(lg::LifeGrid, buffer, I, J, sparseordense)

Update this grid chunk and its associated halos, and return its last cluster.

The last cluster, which is stored before the update, is used as the first cluster in the
buffer for the update of the next chunk down in the `J`th column.
"""
Base.@propagate_inbounds @inline function updatebufferedchunk!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    buffer::Vector{C},
    I::Integer,
    J::Integer,
    ::Type{Sp},
) where {R,C,H,Tall,Wide,Sp<:Sparsity}
    previous = buffer[end-1]
    updategridchunk!(lg, buffer, I, J, Sp)
    updatehalos!(lg, I, J)
    return previous
end



"""
    updatecolumn!(lg::LifeGrid, J, threadid, sparseordense)

Update the `J`th column of `lg.grid` in place.

`threadid` is used to deconflict buffer use between threads; no two instances of
`updatecolumn!` should ever have the same `threadid` simultaneously.

`sparseordense` can be `sparse` or `dense`. The `sparse` algorithm uses a little bit of
overhead (very little for large grids, several % for small ones) to allow it to skip over
fully inactive chunks.
"""
Base.@propagate_inbounds @inline function updatecolumn!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    J::Integer,
    threadid = 1,
) where {R,C,H,Tall,Wide}
    return updatecolumn!(lg, J, threadid, dense)
end

Base.@propagate_inbounds @inline function updatecolumn!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    J::Integer,
    threadid,
    ::Type{dense},
) where {R,C,H,Tall,Wide}
    currentbuffer = lg.buffers[2*threadid-1]
    nextbuffer = lg.buffers[2*threadid]
    previous = zero(C)
    Irange = axes(lg.halos.currentleft, 1)

    updatebuffers!(currentbuffer, lg, first(Irange), J)

    for I in Irange
        currentbuffer[begin] = previous
        if I == last(Irange)
            currentbuffer[end] = zero(C)
        else
            updatebuffers!(nextbuffer, lg, I + 1, J)
            currentbuffer[end] = nextbuffer[2]
        end

        previous = updatebufferedchunk!(lg, currentbuffer, I, J, dense)
        currentbuffer, nextbuffer = nextbuffer, currentbuffer
    end

    return nothing
end

Base.@propagate_inbounds @inline function updatecolumn!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    J::Integer,
    threadid,
    ::Type{sparse},
) where {R,C,H,Tall,Wide}
    buffer = lg.buffers[2*threadid-1]
    changed = lg.changed.current
    previous = zero(C)
    haveprevious = false

    for I in axes(lg.halos.currentleft, 1)
        if !chunkactive(changed, I, J)
            lg.changed.next[I + 1, J + 1] = zero(UInt16)
            haveprevious = false
            continue
        end

        updatebuffers!(buffer, lg, I, J)
        buffer[begin] = haveprevious ? previous : haloedcluster(lg, I - 1, J, nbits(H))
        buffer[end] = haloedcluster(lg, I + 1, J, 1)

        previous = updatebufferedchunk!(lg, buffer, I, J, sparse)
        haveprevious = true
    end

    return nothing
end
