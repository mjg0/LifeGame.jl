export step!, Serial, Parallel



"""
    updatedchunkhalos(chunk::AbstractVector, ::Type{H})

Return the outgoing left and right halo cells of `chunk` as packed bits in two `H`s.

The second bit in the first element of `chunk` goes into the top bit of the left halo, the
second bit of the next element into the second bit of the left halo, etc. Likewise, the
second from the last bit in the first element of `chunk` goes into the top bit of the right
halo, etc. `chunk` must thus be of length `8*sizeof(H)`.
"""
Base.@propagate_inbounds @inline function updatedchunkhalos(
    chunk::AbstractVector{C},
    ::Type{H},
) where {C,H}
    @boundscheck if length(chunk) != nbits(H)
        throw(ArgumentError("chunk must be of length 8*sizeof($H)"))
    end

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

    return lhalo, rhalo
end



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
    stepraw!(lg::LifeGrid)

Advance `lg` by one generation in place and return it.

There are 3 sets of matrices involved in updating `lg`:

1. `lg.data`:

Each column is split into chunks of length equal to the number of bits in the halos.
Boundaries are special cases, but the rest of the update descends through each column thus:

1. The following chunk's clusters are placed into a buffer with their halo bits updated
1. The current chunk's clusters are updated from the previous iteration's buffer
1. The newly updated chunk's halo bits are extracted for future iterations
"""



function swapbuffers(buffers)
    buffers.current, buffers.next = buffers.next, buffers.current
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
    I1, I2 = 1, gridheight(lg)
    J1, J2 = 1, gridheight(lg) # compensate for ghost columns at edges

    @inbounds @batch_if Par for J = J1:J2
        buffer = lg.buffers[Threads.threadid()]

        updatebuffers!(buffer, lg, 1, J)
        # grid[1,J] -> buffer[1,J]

        for I = I1:(I2-1)
            updatebuffers!(buffer, lg, I + 1, J)
            # grid[I,J], buffer[I,J] -> buffer[I,J], next[I+1,J]

            updategridchunk!(lg, buffer, I, J)
            # buffer[I,J] -> grid[I,J]

            updatehalos!(lg, I, J)
            # grid[I,J] -> halos[I,J]
        end

        swapbuffers(buffer)

        #buffer.current, buffer.next = buffer.next, buffer.current
        updategridchunk!(lg, buffer, I2, J)
        # buffer[I,J] -> grid[I,J]

        updatehalos!(lg, I2, J)
    end

    # Swap halos
    lg.halos.currentleft, lg.halos.nextleft = lg.halos.nextleft, lg.halos.currentleft
    lg.halos.currentright, lg.halos.nextright = lg.halos.nextright, lg.halos.currentright

    return lg
end

function step!(lg::LifeGrid, threadmode::Symbol)
    if threadmode === :serial
        return step!(lg, Serial)
    elseif threadmode === :parallel
        return step!(lg, Parallel)
    else
        throw(ArgumentError("threadmode must be either :serial or :parallel"))
    end
end

step!(lg::LifeGrid) = step!(lg, size(lg.grid, 2) > 3 ? Parallel : Serial)



gridheight(lg::LifeGrid) = size(lg.halos.currentleft, 1)

gridwidth(lg::LifeGrid) = size(lg.halos.currentleft, 2) - 2



Base.@propagate_inbounds @inline function gridchunk(
    lg::LifeGrid{R,C,H,Tall,Wide},
    I,
    J,
) where {R,C,H,Tall,Wide}
    i = (I - 1) * nbits(H) + 2
    return view(lg.grid, i:(i+nbits(H)-1), J + 1)
end



"""
Begin state: buffer.next has the previous iteration's results
End state: buffer.current has this iteration's result

Reads from lg.grid, lg.halos, and buffer.next, writes to buffer.current and buffer.next; swaps buffer.{current,next}
"""
Base.@propagate_inbounds @inline function updatebuffers!(
    buffer::Buffer{C},
    lg::LifeGrid{R,C,H,Tall,Wide},
    I::Integer,
    J::Integer,
) where {R,C,H,Tall,Wide}
    # Convenience variables
    chunk = gridchunk(lg, I, J)
    cbuf = buffer.current
    nbuf = buffer.next
    lhalo = lg.halos.currentright[I, J]
    rhalo = lg.halos.currentleft[I, J+2]

    # Before overwriting it, store the last real element of cbuf
    nbuf[begin] = ifelse(I == 1, zero(C), cbuf[end-1])

    #if Wide
        # Apply lhalo and rhalo to chunk, storing the results in nbuf
        @simd for k = 1:nbits(H)
            centermask = ~(lowbit(C) | highbit(C)) # all cells but the outermost two on
            lbit = C(lhalo >> (nbits(H) - k) & one(H)) << (nbits(C) - 1)
            rbit = C(rhalo >> (nbits(H) - k) & one(H))

            cbuf[k+1] = (chunk[k] & centermask) | lbit | rbit
        end
    #else
        # If there are no halos that require updating, a copy is all that's needed
        cbuf[(begin+1):(end-1)] .= chunk
    #end

    # Get the first real element of nbuf in preparation for the next iteration
    nbuf[end] = ifelse(I == gridheight(lg), zero(C), cbuf[begin+1])

    buffer.current, buffer.next = buffer.next, buffer.current

    return nothing
end



"""
Reads from buffers.current, writes to lg.grid
"""
Base.@propagate_inbounds @inline function updategridchunk!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    buffer::Buffer{C},
    I::Integer,
    J::Integer,
) where {R,C,H,Tall,Wide}
    chunk = gridchunk(lg, I, J)
    cbuf = buffer.current

    # Update each cluster in chunk
    @simd for i = 1:nbits(H)
        chunk[i] = updatedcluster(cbuf[i], cbuf[i+1], cbuf[i+2], R)
    end

    # Zero out trailing columns
    if J == gridwidth(lg)-1
        # How many bits to zero out
        n = mod(-size(lg, 2), nbits(C) - 2) + 1

        # Roll off n bits for each cluster
        @simd for i = 1:nbits(H)
            chunk[i] = (chunk[i] >> n) << n
        end
    end

    # Zero out trailing rows
    if I == gridheight(lg)
        lg.grid[lg.height+2, J] = zero(C)
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

    lg.halos.nextleft[I, J+1] = lhalo
    lg.halos.nextright[I, J+1] = rhalo

    return nothing
end

updatehalos!(::LifeGrid{R,C,H,Tall,false}, args...) where {R,C,H,Tall} = nothing



#Base.@propagate_inbounds @inline function updatehalos!(
#    out::AbstractVector{T},
#    in::AbstractVector{T},
#    lhalo::H,
#    rhalo::H,
#) where {T,H}
#    centermask = ~(lowbit(T) | highbit(T))
#
#    @simd for k = 1:nbits(H)
#        lbit = T(lhalo >> (nbits(H) - k) & one(H)) << (nbits(T) - 1)
#        rbit = T(rhalo >> (nbits(H) - k) & one(H))
#
#        out[k+1] = (in[k] & centermask) | lbit | rbit
#    end
#end
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#"""
#
#`out` is 2 longer than `in` and doesn't get its first or last cells updated
#"""
#Base.@propagate_inbounds @inline function updatehalos!(
#    out::AbstractVector{T},
#    in::AbstractVector{T},
#    lhalo::H,
#    rhalo::H,
#) where {T,H}
#    centermask = ~(lowbit(T) | highbit(T))
#
#    @simd for k = 1:nbits(H)
#        lbit = T(lhalo >> (nbits(H) - k) & one(H)) << (nbits(T) - 1)
#        rbit = T(rhalo >> (nbits(H) - k) & one(H))
#
#        out[k+1] = (in[k] & centermask) | lbit | rbit
#    end
#end
#
#
#
#
#
#Base.@propagate_inbounds @inline function updatehalos!(lg::LifeGame{R,C,H,Wide,Tall}, I::Integer, J::Integer) where {R,C,H,Wide,Tall}
#    centermask = ~(lowbit(T) | highbit(T))
#
#    chunk = gridchunk(lg, I, J)
#    buffer = lg.buffers.current
#
#    if Wide
#        lhalo = lg.halos.currentright[I,J]
#        rhalo = lg.halos.currentleft[I,J]
#
#        @simd for k = 1:nbits(H)
#            lbit = T(lhalo >> (nbits(H) - k) & one(H)) << (nbits(T) - 1)
#            rbit = T(rhalo >> (nbits(H) - k) & one(H))
#
#            buffer[k+1] = (chunk[k] & centermask) | lbit | rbit
#        end
#    else
#        buffer[begin+1:end-1]
#    end
#end
#
#
#
#
#"""
#
#`in` should be two elements longer than `out`
#"""
#Base.@propagate_inbounds @inline function updategridchunk!(
#    out::AbstractVector,
#    in::AbstractVector,
#    rule::R,
#    ::Type{H},
#) where {R,H}
#    @simd for i = 1:(8*sizeof(H))
#        out[i] = updatedcluster(in[i], in[i+1], in[i+2], rule)
#    end
#end
#
#
#
#
#
#
#Base.@propagate_inbounds @inline function zerotrailing!(
#    cells::AbstractVector,
#    shift,
#    ::Type{H},
#) where {H}
#    @simd for i = 1:(8*sizeof(H))
#        cells[i] = (cells[i] >> shift) << shift
#    end
#end
#
