export step!



"""
    updatedhalos(left::Integer, current::Integer, right::Integer)

Return `current` with its leftmost and rightmost bit updated given `left` and `right`.

Flips the first bit of `current` to the value of the second-to-last bit of `left` (`left`'s
last active bit), and the last bit of `current` to the value of the second bit of `right`
(`right`'s first active bit).
"""
function updatedhalos(left::T, current::T, right::T) where T
    N = 8*sizeof(T)
    firstbit = one(T) << (N-1)
    lastbit = one(T)

    lefthalo = (left & (lastbit    << 1)) << (N-2)
    righthalo = (right & (firstbit >> 1)) >> (N-2)
    return (current & ~(lastbit | firstbit)) | lefthalo | righthalo
end

Base.@propagate_inbounds @inline function updatedhalos(cluster::T, lefthalo::H, righthalo::H, ::Val{K}) where {T,H,K}
    Hbits = 8*sizeof(H)
    Tbits = 8*sizeof(T)

    centermask = (typemax(T) << 2) >> 1
    lbit = T(lefthalo  >> (Hbits-K) & one(H)) << (Tbits-1)
    rbit = T(righthalo >> (Hbits-K) & one(H))

    return (cluster & centermask) | lbit | rbit
end

Base.@propagate_inbounds @inline function updatedhalos(cluster::T, lefthalo::H, righthalo::H, k::Integer) where {T,H,K}
    Hbits = 8*sizeof(H)
    Tbits = 8*sizeof(T)

    centermask = (typemax(T) << 2) >> 1
    lbit = T(lefthalo  >> (Hbits-k) & one(H)) << (Tbits-1)
    rbit = T(righthalo >> (Hbits-k) & one(H))

    return (cluster & centermask) | lbit | rbit
end

#@inline function updatebdchunkhalos(block::NTuple{N,T}, ::Type{H}) where {N,T<:Unsigned,H<:Unsigned}
Base.@propagate_inbounds @inline function updatedchunkhalos(lg::LifeGrid{R,C,H}, i, j) where {R,C,H}
    lhalo = zero(H)
    rhalo = zero(H)

    lshift = 8*sizeof(C)-2
    rshift = 1

    N = 8*sizeof(H)
    for k in 1:N
        lhalo |= H((lg.grid[i+k-1, j] >> lshift) & one(C)) << (N-k)
        rhalo |= H((lg.grid[i+k-1, j] >> rshift) & one(C)) << (N-k)
    end

    return lhalo, rhalo
end
#Base.@propagate_inbounds @inline @generated function updatedchunkhalos(grid::AbstractMatrix{T}, i, j, ::Type{H}) where {T,H}
#    N = 8*sizeof(H)
#
#    lshift = 8*sizeof(T)-2
#    rshift = 1
#
#    updates = Expr(:block, Iterators.flatten((
#        (
#            :(lhalo |= H((grid[i+$(k-1), j] >> $lshift) & one(T)) << $(N-k)),
#            :(rhalo |= H((grid[i+$(k-1), j] >> $rshift) & one(T)) << $(N-k)),
#        ) for k in 1:N
#    ))...)
#
#    return quote
#        lhalo = zero(H)
#        rhalo = zero(H)
#
#        $(updates)
#
#        return lhalo, rhalo
#    end
#end



"""
    step!(lg::LifeGrid; chunklength=$DEFAULT_CHUNK_SIZE, parallel=size(lg, 1)>1024)

Update `lg` one generation according to the [`rule`](@ref) associated with it.

All cells outside of the grid boundary are fixed at zero.

`step!` runs using all available threads by default.

`chunklength` determines the size of a chunk of data that `step!` works on before proceeding
to the next chunk. $DEFAULT_CHUNK_SIZE is chosen as the default since it strikes a good
balance: it leads to chunks large enough that work isn't interrupted too often, and small
enough to fit in the L1 cache of most machines. The height of `lg` must be at least
`chunksize*Threads.nthreads()` for all threads to be fully engaged.

`parallel` determines whether `step!` will run with multiple threads. It is `true` by
default if `lg`'s height exceeds 1024, and `false` otherwise. This is a reasonable default
on most machines, but it's worth experimenting with.

A generic algorithm for updating each cluster in the grid is used by default. The compiler
does a decent job of optimizing for most rules, but hand-tuning the cluster update function
can improve performance by 10% for some rules. See the extended help for
[`LifeGame.updatedcluster`](@ref) for instructions on specializing the cluster update.
Specializations are provided for commonly used rules (`B3/S23`, `B36/S23`, and `B2/s`).
"""
function step!(lg::LifeGrid; chunklength=DEFAULT_CHUNK_SIZE, parallel=size(lg, 1)>1024)
    if parallel
        return stepraw!(lg, chunklength)
    else
        disable_polyester_threads() do
            return stepraw!(lg, chunklength)
        end
    end
end



#"""
#    stepraw!(lg::LifeGrid, chunklength)
#
#The back end of [`step!`](@ref); updates and returns `lg`.
#
## Extended help
#
#`lg` is updated one column of clusters at a time--that is, a column of 64-bit unsigned
#integers, each representing a 62-cell portion of a row (see [`updatedcluster`](@ref)).
#
#There are 3 arrays involved in the update:
#
#1. The grid itself, a `Matrix` of `UInt64` clusters.
#2. The left buffer, a `Vector` of `UInt64`s with as many elements as the grid has rows.
#3. The middle buffer, identical in type and size to the left buffer.
#
#These are carefully updated in an order that eliminates inter-iteration dependencies when
#updating a column. This allows the compiler to emit vectorized instructions and makes it
#safe to operate on a single column with multiple threads. The use of only two buffer rows,
#rather than an entire intermediate grid, keeps the memory footprint small and makes better
#use of the cache.
#
#Here is how one cluster on the interior of the grid is updated, following the actual order of
#operations used in this function:
#
#1. The cluster is updated to according to the life rules using the corresponding cells
#   above, at the same place as, and below the cluster in the left buffer using
#   `updatedcluster`:\n
#   `grid[i,j] = updatedcluster(leftbuffer[i-1], leftbuffer[i], leftbuffer[i+1])`.
#1. The corresponding cluster in the middle buffer has its halos updated; the corresponding
#   cluster in the left buffer is used as the left value, the cluster to the right of the
#   grid cluster in question as the middle value, and the cluster two to the right of the
#   grid cluster in question as the right value:\n
#   `rightbuffer[i] = updatedhalos(leftbuffer[i], grid[i,j+1], grid[i,j+2])`.
#1. The buffers are switched in preparation for the next iteration.
#
#Notice that no cell that is written to on this iteration is also read from; this allows for
#vectorization and multithreading. Cache efficiency is increased by doing these updates in
#chunks: taking a portion of a column, updating it fully, and only then moving on to the next
#chunk.
#
#Columns at the boundaries are special cases, but the order of operations is the same. The
#computations are aided by having padding columns to the left and right of the first and last
#active grids.
#"""
Base.@propagate_inbounds @inline @generated function halotuple(grid::AbstractMatrix{T}, lhalo::H, rhalo::H, i, j) where {T,H}
    return Expr(:tuple, (
        :(updatedhalos(grid[i+$(k-1),j], lhalo, rhalo, Val($k)))
        for k in 1:8*sizeof(H)
    )...)
end

Base.@propagate_inbounds @inline function halotuple!(buffer::AbstractVector{T}, grid::AbstractMatrix{T}, lhalo::H, rhalo::H, i, j) where {T,H}
    N = 8*sizeof(H)
    @simd for k in 1:N
        buffer[k+1] = updatedhalos(grid[i+k-1,j], lhalo, rhalo, k)
    end
end



#Base.@propagate_inbounds @inline @generated function updategridchunk!(grid::AbstractMatrix{T}, i, j, above0, cur::NTuple{N,T}, belowB, rule) where {N,T}
#    return Expr(:block, (
#        :(grid[i,j] = updatedcluster(above0, cur[1], cur[2], rule)),
#        (
#            :(grid[i+$(k-1),j] = updatedcluster(cur[$(k-1)], cur[$k], cur[$(k+1)], rule))
#            for k in 2:N-1
#        )...,
#        :(grid[i+$(N-1),j] = updatedcluster(cur[$(N-1)], cur[$N], belowB, rule)),
#        :(return nothing)
#    )...)
#end

Base.@propagate_inbounds @inline function updategridchunk!(lg::LifeGrid{R,C,H}, i, j, cur::AbstractVector{C}) where {R,C,H}
    N = 8*sizeof(H)
    @simd for k in 1:(N-2)
        lg.grid[i+k-1,j] = updatedcluster(cur[k], cur[k+1], cur[k+2], R)
    end
end

Base.@propagate_inbounds @inline @generated function padtuple(t::NTuple{N,T}, front::T, back::T)::NTuple{N+2,T} where {N,T}
    return Expr(:tuple,
                :(front),
                (:(t[$i]) for i in 1:N)...,
                :(back))
end



function stepraw!(lg::LifeGrid{R,C,H}, dummy) where {R,C,H}
    Hbits = 8*sizeof(H)

    grid = lg.grid

    J1 = firstindex(lg.grid, 2)+1
    J2 = lastindex(lg.grid, 2)-1

    bufflen = Hbits+2

    @inbounds @batch per=thread for j in J1:J2
        current = lg.colbuffers1[Threads.threadid()]
        next = lg.colbuffers2[Threads.threadid()]
        inhalosleft   = view(lg.inhalosright,  :, j-1)
        inhalosright  = view(lg.inhalosleft,   :, j+1)
        outhalosleft  = view(lg.outhalosleft,  :, j)
        outhalosright = view(lg.outhalosright, :, j)

        above = zero(C)
        halotuple!(current, grid, inhalosleft[begin], inhalosright[begin], 2, j)
        #current = halotuple(grid, inhalosleft[begin], inhalosright[begin], 2, j)

        # Update this row
        for I in firstindex(inhalosleft):lastindex(inhalosleft)-1
            i = (I-1)*Hbits+2

            #next = halotuple(grid, inhalosleft[I+1], inhalosright[I+1], i+Hbits, j)
            halotuple!(next, grid, inhalosleft[I+1], inhalosright[I+1], i+Hbits, j)
            current[1] = above
            current[bufflen] = next[2]

            updategridchunk!(lg, i, j, current)
            outhalosleft[I], outhalosright[I] = updatedchunkhalos(lg, i, j)

            above = current[bufflen-1]
            current, next = next, current
        end

        I = lastindex(inhalosleft)
        i = (I-1)*Hbits+2

        current[1] = above
        current[bufflen] = zero(C)
        updategridchunk!(lg, i, j, current)
        outhalosleft[I], outhalosright[I] = updatedchunkhalos(lg, i, j)

        # Zero the last+1 cell
        grid[lg.height+1,j] = zero(C)
    end

    return lg
end