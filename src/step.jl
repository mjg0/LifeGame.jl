export step!



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
function step!(lg::LifeGrid{R,C,H}, ::Type{Par}) where {R,C,H,Par<:ThreadingMode}
    grid = lg.grid

    J1 = firstindex(lg.grid, 2) + 1
    J2 = lastindex(lg.grid, 2) - 1

    bufflen = nbits(H) + 2

    endshift = mod(-size(lg, 2), nbits(C) - 2) + 1

    @inbounds @batch_if Par for J = J1:J2
        current = lg.buffers.a[Threads.threadid()]
        next = lg.buffers.b[Threads.threadid()]

        inhalosleft = view(lg.halos.currentright, :, J - 1)
        inhalosright = view(lg.halos.currentleft, :, J + 1)
        outhalosleft = view(lg.halos.nextleft, :, J)
        outhalosright = view(lg.halos.nextright, :, J)

        above = zero(C)
        updatehalos!(current, gridchunk(lg, 1, J), inhalosleft[1], inhalosright[1])

        # Update this row
        for I = firstindex(inhalosleft):(lastindex(inhalosleft)-1)
            updatehalos!(next, gridchunk(lg, I + 1, J), inhalosleft[I+1], inhalosright[I+1])

            current[1] = above
            current[bufflen] = next[2]

            updategridchunk!(gridchunk(lg, I, J), current, R, H)

            # Zero trailing cells
            if J == J2
                zerotrailing!(gridchunk(lg, I, J), endshift, H)
            end

            outhalosleft[I], outhalosright[I] = updatedchunkhalos(gridchunk(lg, I, J), H)

            above = current[bufflen-1]

            current, next = next, current
        end

        I = lastindex(inhalosleft)

        current[1] = above
        current[bufflen] = zero(C)

        updategridchunk!(gridchunk(lg, I, J), current, R, H)

        if J == J2
            zerotrailing!(gridchunk(lg, I, J), endshift, H)
        end

        # Zero the last+1 cell
        grid[lg.height+2, J] = zero(C)

        outhalosleft[I], outhalosright[I] = updatedchunkhalos(gridchunk(lg, I, J), H)
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

step!(lg::LifeGrid) = step!(lg, size(lg.grid, 2)>3 ? Parallel : Serial)



"""

`out` is 2 longer than `in` and doesn't get its first or last cells updated
"""
Base.@propagate_inbounds @inline function updatehalos!(
    out::AbstractVector{T},
    in::AbstractVector{T},
    lhalo::H,
    rhalo::H,
) where {T,H}
    centermask = ~(lowbit(T) | highbit(T))

    @simd for i = 1:nbits(H)
        lbit = T(lhalo >> (nbits(H) - i) & one(H)) << (nbits(T) - 1)
        rbit = T(rhalo >> (nbits(H) - i) & one(H))

        out[i+1] = (in[i] & centermask) | lbit | rbit
    end
end




"""

`in` should be two elements longer than `out`
"""
Base.@propagate_inbounds @inline function updategridchunk!(
    out::AbstractVector,
    in::AbstractVector,
    rule::R,
    ::Type{H},
) where {R,H}
    @simd for i = 1:(8*sizeof(H))
        out[i] = updatedcluster(in[i], in[i+1], in[i+2], rule)
    end
end



Base.@propagate_inbounds @inline function gridchunk(lg::LifeGrid{R,C,H}, I, J) where {R,C,H}
    i = (I - 1) * nbits(H) + 2
    return view(lg.grid, i:(i+nbits(H)-1), J)
end



Base.@propagate_inbounds @inline function zerotrailing!(
    cells::AbstractVector,
    shift,
    ::Type{H},
) where {H}
    @simd for i = 1:(8*sizeof(H))
        cells[i] = (cells[i] >> shift) << shift
    end
end
