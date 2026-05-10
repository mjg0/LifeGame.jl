export step!, serial, parallel



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

    # Convenience variables
    chunk = gridchunk(lg, I, J)

    # Incoming halos are the ones from adjacent columns
    lhalo = J == 1 ? zero(C) : lg.halos.currentright[I, J-1]
    rhalo = J == size(lg.grid, 2) ? zero(C) : lg.halos.currentleft[I, J+1]

    # Apply lhalo and rhalo to chunk, storing the results in buffer
    @simd for k = 1:nbits(H)
        # Mask and bits to update this cluster
        centermask = ~(lowbit(C) | highbit(C)) # all cells but the outermost two on
        lbit = C(lhalo >> (nbits(H) - k) & one(H)) << (nbits(C) - 1)
        rbit = C(rhalo >> (nbits(H) - k) & one(H))

        # Put the updated cluster in buffer
        if Wide
            buffer[k+1] = (chunk[k] & centermask) | lbit | rbit
        else
            buffer[k+1] = chunk[k] & centermask
        end
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
    chunk = gridchunk(lg, I, J)
    cbuf = buffer

    # Update each cluster in chunk
    @simd for i = 1:nbits(H)
        chunk[i] = updatedcluster(cbuf[i], cbuf[i+1], cbuf[i+2], R)
    end

    # Zero out trailing columns
    if J == size(lg.grid, 2)
        # How many bits to zero out
        n = mod(-size(lg, 2), nbits(C) - 2) + 1

        # Roll off n bits for each cluster
        @simd for k = 1:nbits(H)
            chunk[k] = (chunk[k] >> n) << n
        end
    end

    # Zero out trailing rows
    if I == size(lg.halos.currentleft, 1) && size(lg.grid, 1) > lg.height
        lg.grid[lg.height+1, J] = zero(C)
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
    updatecolumn!(lg::LifeGrid, J, tid=1)

Update the `J`th column of `lg.grid` in place.

`tid` is used to deconflict buffer use between threads; no two instances of `updatecolumn!`
should ever have the same `tid` simultaneously.

# Extended help

Several arrays are involved in the update of a column of clusters:

- `lg.grid`: the main piece of `lg`, a matrix of clusters, unsigned integers whose packed
  bits each represent a single cell. Each of these clusters represents a portion of a row.
  It's logically divided into chunks, with each chunk being of length equal to the number of
  bits in `lg`'s halos. For large grids, with both cluster and halo types set to `UInt64`, a
  chunk is 64 `UInt64`s long and represents a 64×62 cell section of the logical grid.
- `lg.halos`: a struct containing 4 matrices, each with as many elements as `lg` has chunks.
  There are left and right "incoming" and left and right "outgoing" halos. At the end of
  [`step!`](@ref), the incoming and outgoing halos are swapped.
- `lg.buffers`: several vectors, each the size of a chunk in `lg.grid` plus padding cells at
  the begining and end. There are enough of these buffers for each thread to use two.

The update for the `(I, J)`th chunk of `lg` proceeds as follows, ignoring minor nuances and
cases at the boundary:

1. The buffer representing the next chunk takes on the values from the `(I, J+1)`th chunk of
   `lg.grid`, with the clusters' halo bits updated from the incoming left and right halos.
2. The `(I, J)`th chunk of `lg.grid` has its clusters fully updated, reading from the buffer
   representing the current chunk, which had its halo bits correctly set last iteration.
3. The outgoing halos are computed from the bits in the newly updated `(I, J)`th chunk of
   `lg.grid` and stored in anticipation of the next call to `step!`.
"""
Base.@propagate_inbounds @inline function updatecolumn!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    J::Integer,
    threadid = 1,
) where {R,C,H,Tall,Wide}
    # Select two unique buffers
    currentbuffer = lg.buffers[2*threadid-1]
    nextbuffer = lg.buffers[2*threadid]

    # Initialize the first buffer and the carry variable
    updatebuffers!(currentbuffer, lg, 1, J)
    previous = zero(C)

    # Loop over all the the last chunk in this column
    Irange = axes(lg.halos.currentleft, 1)
    for I = first(Irange):(last(Irange)-1)
        # Skip one chunk ahead to update the next buffer
        currentbuffer[begin] = previous
        updatebuffers!(nextbuffer, lg, I + 1, J)
        currentbuffer[end] = nextbuffer[2]

        # Update this chunk from currentbuffer
        updategridchunk!(lg, currentbuffer, I, J)

        # Calculate halos based off of the update
        updatehalos!(lg, I, J)

        # Switch buffers
        previous = currentbuffer[end-1]
        currentbuffer, nextbuffer = nextbuffer, currentbuffer
    end

    # Update the first and last cells in currentbuffer
    currentbuffer[begin] = previous
    currentbuffer[end] = zero(C)

    # Update the last chunk
    updategridchunk!(lg, currentbuffer, last(Irange), J)

    # Last halos in this column
    updatehalos!(lg, last(Irange), J)

    return nothing
end



# @batch is noticeably expensive for small grids, this allows that to be avoided
abstract type ThreadingMode end
struct serial <: ThreadingMode end
struct parallel <: ThreadingMode end



"""
    step!(lg::LifeGrid, threadmode)

Update `lg` one generation according to the [`rule`](@ref) associated with it.

All cells outside of the grid boundary are fixed at zero.

`threadmode` can be `serial` or `parallel`. By default, it's `parallel` if the grid backing
`lg` has multiple columns.
"""
@inline function step!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    ::Type{Par},
) where {R,C,H,Tall,Wide,Par<:ThreadingMode}
    # Only parallelize if there's a need
    if Wide && Par === parallel
        # Manual thread bounds so buffers don't get mixed
        nthreads = min(Threads.nthreads(), size(lg.grid, 2))
        colsperthread = cld(size(lg.grid, 2), nthreads)

        # Outer loop over threads, inner over columns
        @inbounds @batch for tid = 1:nthreads
            J1 = colsperthread * (tid - 1) + 1
            J2 = min(colsperthread * tid, size(lg.grid, 2))
            for J = J1:J2
                updatecolumn!(lg, J, tid)
            end
        end
    else
        @inbounds for J in axes(lg.grid, 2)
            updatecolumn!(lg, J)
        end
    end

    # Swap halos
    lg.halos.currentleft, lg.halos.nextleft = lg.halos.nextleft, lg.halos.currentleft
    lg.halos.currentright, lg.halos.nextright = lg.halos.nextright, lg.halos.currentright

    return lg
end

step!(lg::LifeGrid) = step!(lg, parallel)