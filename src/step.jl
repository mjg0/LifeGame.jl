export step!



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



# disable_polyester_threads is noticeably expensive for small grids, @batch_if allays that
abstract type ThreadingMode end
struct Serial <: ThreadingMode end
struct Parallel <: ThreadingMode end
macro batch_if(mode, loop)
    return esc(quote
        if $mode === Parallel
            @batch $loop
        elseif $mode === Serial
            $loop
        else
            throw(ArgumentError("expected Serial or Parallel"))
        end
    end)
end



"""
    step!(lg::LifeGrid)
    step!(lg::LifeGrid, threadmode::Symbol)

Update `lg` one generation according to the [`rule`](@ref) associated with it.

All cells outside of the grid boundary are fixed at zero.

`threadmode` can be `:serial` or `:parallel`. By default, it's `:parallel` only if the grid
backing `lg` has multiple columns.
"""
function step!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    ::Type{Par},
) where {R,C,H,Tall,Wide,Par<:ThreadingMode}
    # Get iteration bounds
    I1, I2 = 1, size(lg.halos.currentleft, 1)

    # @batch doesn't like "end"
    bufflen = nbits(H) + 2

    # Manual thread bounds so buffers don't get mixed
    nthreads = min(Threads.nthreads(), size(lg.grid, 2))
    colsperthread = cld(size(lg.grid, 2), nthreads)

    @inbounds @batch_if Par for tid in 1:nthreads
        J1 = colsperthread * (tid - 1) + 1
        J2 = min(colsperthread * tid, size(lg.grid, 2))

        for J in J1:J2
            # Select two unique buffers
            tid = Threads.threadid()
            currentbuffer = lg.buffers[2*tid-1]
            nextbuffer = lg.buffers[2*tid]

            # Initialize the first buffer
            updatebuffers!(currentbuffer, lg, 1, J)

            previous = zero(C)
            for I = I1:(I2-1)
                # Skip one chunk ahead to update the next buffer
                currentbuffer[1] = previous
                updatebuffers!(nextbuffer, lg, I + 1, J)
                currentbuffer[bufflen] = nextbuffer[2]

                # Update this chunk from currentbuffer
                updategridchunk!(lg, currentbuffer, I, J)

                # Calculate halos based off of the update
                updatehalos!(lg, I, J)

                # Switch buffers and 
                previous = currentbuffer[bufflen-1]
                currentbuffer, nextbuffer = nextbuffer, currentbuffer
            end

            # Update the first and last cells in currentbuffer
            currentbuffer[1] = previous
            currentbuffer[bufflen] = zero(C)

            # Update the last chunk
            updategridchunk!(lg, currentbuffer, I2, J)

            # Last halos in this column
            updatehalos!(lg, I2, J)
        end
    end

    # Swap halos
    lg.halos.currentleft, lg.halos.nextleft = lg.halos.nextleft, lg.halos.currentleft
    lg.halos.currentright, lg.halos.nextright = lg.halos.nextright, lg.halos.currentright

    return lg
end

function step!(lg::LifeGrid{R,C,H,Tall,Wide}, threadmode::Symbol) where {R,C,H,Tall,Wide}
    if !Wide || threadmode === :serial
        return step!(lg, Serial)
    elseif threadmode === :parallel
        return step!(lg, Parallel)
    else
        throw(ArgumentError("threadmode must be either :serial or :parallel"))
    end
end

step!(lg::LifeGrid) = step!(lg, :parallel)