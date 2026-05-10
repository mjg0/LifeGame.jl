export step!



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

`threadmode` can be `:serial` or `:parallel`. By default, it's set to true if the grid
backing `lg` has multiple columns.

A generic algorithm for updating each cluster in the grid is used by default. The compiler
does a decent job of optimizing for most rules, but hand-tuning the cluster update function
can improve performance by 10% for some rules. See the extended help for
[`LifeGame.updatedcluster`](@ref) for instructions on specializing the cluster update.
Specializations are provided for commonly used rules (`B3/S23`, `B36/S23`, and `B2/s`).
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
        J = J1:J2
        for J in J1:J2
            # Select two unique buffers
            tid = Threads.threadid()
            currentbuffer = lg.buffers[2*tid-1]
            nextbuffer = lg.buffers[2*tid]

            # Initialize the first buffer
            updatebuffers!(currentbuffer, lg, 1, J)

            previous = zero(C)
            for I = I1:(I2-1)
                # Skip one step ahead to the next buffer
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



Base.@propagate_inbounds @inline function gridchunk(
    lg::LifeGrid{R,C,H,Tall,Wide},
    I,
    J,
) where {R,C,H,Tall,Wide}
    i = (I - 1) * nbits(H) + 1
    return view(lg.grid, i:(i+nbits(H)-1), J)
end



"""
Begin state: buffer.next has the previous iteration's result

End state: buffer.current has this iteration's result

Reads from lg.grid, lg.halos, and buffer.next, writes to buffer.current and buffer.next; swaps buffer.{current,next}
"""
Base.@propagate_inbounds @inline function updatebuffers!(
    buffer::Vector{C},
    lg::LifeGrid{R,C,H,Tall,Wide},
    I::Integer,
    J::Integer,
) where {R,C,H,Tall,Wide}
    # Convenience variables
    chunk = gridchunk(lg, I, J)

    # Incoming halos are the ones from adjacent columns
    lhalo = J == 1 ? zero(C) : lg.halos.currentright[I, J-1]
    rhalo = J == size(lg.grid, 2) ? zero(C) : lg.halos.currentleft[I, J+1]

    # Apply lhalo and rhalo to chunk, storing the results in buffer
    @simd for k = 1:nbits(H)
        centermask = ~(lowbit(C) | highbit(C)) # all cells but the outermost two on
        lbit = C(lhalo >> (nbits(H) - k) & one(H)) << (nbits(C) - 1)
        rbit = C(rhalo >> (nbits(H) - k) & one(H))

        if Wide
            buffer[k+1] = (chunk[k] & centermask) | lbit | rbit
        else
            buffer[k+1] = chunk[k] & centermask
        end
    end

    return nothing
end



"""
Reads from buffers.current, writes to lg.grid
"""
Base.@propagate_inbounds @inline function updategridchunk!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    buffer, #::Buffer{C},
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
Reads from lg.grid and writes to lg.halos
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