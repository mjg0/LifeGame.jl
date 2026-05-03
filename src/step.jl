export step!



"""
    updatedhalos(left::Integer, current::Integer, right::Integer)

Return `current` with its leftmost and rightmost bit updated given `left` and `right`.

Flips the first bit of `current` to the value of the second-to-last bit of `left` (`left`'s
last active bit), and the last bit of `current` to the value of the second bit of `right`
(`right`'s first active bit).
"""
function updatedhalos(left, current, right)
    lefthalo = (left & (LAST_BIT << 1)) << CELLS_PER_CLUSTER
    righthalo = (right & (FIRST_BIT  >> 1)) >> CELLS_PER_CLUSTER
    return (current & ~(LAST_BIT | FIRST_BIT)) | lefthalo | righthalo
end

@inline function updatedhalos(cluster::T, lefthalo::H, righthalo::H, ::Val{K}) where {T,H,K}
    N = 8*sizeof(H)
    centermask = (typemax(T) << 2) >> 1
    lbit = T(lefthalo  >> (N-K) & one(H)) << (8*sizeof(T)-1)
    rbit = T(righthalo >> (N-K) & one(H))
    return (cluster & centermask) | lbit | rbit
end

#@inline function updatebdchunkhalos(block::NTuple{N,T}, ::Type{H}) where {N,T<:Unsigned,H<:Unsigned}
@inline function updatedchunkhalos(grid::AbstractMatrix{T}, i, j, N, ::Type{H}) where {T,H}
    lhalo = zero(H)
    rhalo = zero(H)

    lshift = 8*sizeof(T)-2
    rshift = 1

    @inbounds for k in 1:N
        lhalo |= H((grid[i+k-1,j] >> lshift) & one(T)) << (N-k)
        rhalo |= H((grid[i+k-1,j] >> rshift) & one(T)) << (N-k)
    end

    return lhalo, rhalo
end



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
const _UInt = Unsigned

@inline @inbounds function halotuple(grid, lhalo, rhalo, i, j, ::Val{N}) where N
    return ntuple(k->updatedhalos(grid[i+k-1,j], lhalo, rhalo, Val(k)), Val(N))
end

@generated function updategridchunk!(grid, i, j, above0, cur::NTuple{B,T}, belowB, rule) where {T,B}
    updates = Expr(:block, (
        :(grid[i,j] = updatedcluster(above0, cur[1], cur[2], rule)),
        (
            :(grid[i+$(k-1),j] = updatedcluster(cur[$(k-1)], cur[$k], cur[$(k+1)], rule)
            for k in 2:B-1)
        )...,
        :(grid[i+$(B-1),j] = updatedcluster(cur[$(B-1)], cur[$B], belowB, rule))
    )...)

    return quote
        @inbounds begin
            $updates
        end
        return nothing
    end
end

@generated function updategridendchunk!(grid, i, j, above0, cur::NTuple{M,T}, belowM, rule) where {T,M}
    stores = Expr(:block)

    stores = if M == 1
        Expr(:block, :(grid[i,j] = updatedcluster(above0, cur[i], belowM, rule)))
    else
        Expr(:block, (
            :(grid[i, j] = updatedcluster(above0, cur[1], cur[2], rule)),
            (
                :(grid[i+$(k-1),j] = updatedcluster(cur[$(k-1)], cur[$k], cur[$(k+1)], rule))
                for k in 2:M-1
            )...,
            :(grid[i+$(M-1),j] = updatedcluster(cur[$(M-1)], cur[$M], belowM, rule))
        )...)
    end

    return quote
        @inbounds begin
            $stores
        end
        return nothing
    end
end

function stepraw!(lg::LifeGrid{R}, dummy) where {R}
    grid = lg.grid
    inL  = lg.inhalosleft
    inR  = lg.inhalosright

    H = eltype(inL)
    T = eltype(grid)
    B = 8sizeof(H)

    ilo = firstindex(grid, 1)+1
    ihi = lastindex(grid, 1)-1
    jlo = firstindex(grid, 2)+1
    jhi = lastindex(grid, 2)-1

    ilo > ihi && return lg

    hfirst = firstindex(inL, 1)
    hlast  = lastindex(inL, 1)

    @inbounds @batch for j in jlo:jhi
        i = ilo
        h = hfirst

        above0 = zero(T)

        # Process all chunks except the final chunk.
        while i + B - 1 < ihi
            cur  = halotuple(grid, inL[h],   inR[h],   i,   j, Val(B))
            next = halotuple(grid, inL[h+1], inR[h+1], i+B, j, Val(B))

            updategridchunk!(grid, i, j, above0, cur, next[1], R)
            lg.outhalosleft[h], lg.outhalosright[h] = updatedchunkhalos(grid, i, j, B, H)

            above0 = cur[end]
            i += B
            h += 1
        end

        # Final chunk: may be partial.
        M = ihi - i + 1

        if M == B
            cur = halotuple(grid, inL[h], inR[h], i, j, Val(B))

            # Bottom boundary row. Use real halo entry if present; otherwise zero halos.
            belowB =
                if h < hlast
                    updatedhalos(grid[i+B, j], inL[h+1], inR[h+1], Val(1))
                else
                    updatedhalos(grid[i+B, j], zero(H), zero(H), Val(1))
                end

            updategridchunk!(grid, i, j, above0, cur, belowB, R)
            lg.outhalosleft[h], lg.outhalosright[h] = updatedchunkhalos(grid, i, j, M, H)
        else
            cur = halotuple(grid, inL[h], inR[h], i, j, Val(M))

            # Since this is a partial chunk, the row below is still within the same
            # halo word position, at Val(M+1).
            belowM = updatedhalos(grid[i+M, j], inL[h], inR[h], Val(M+1))

            updategridendchunk!(grid, i, j, above0, cur, belowM, R)
            lg.outhalosleft[h], lg.outhalosright[h] = updatedchunkhalos(grid, i, j, M, H)
        end
    end

    return lg
end