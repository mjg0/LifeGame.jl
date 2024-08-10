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



"""
    step!(lg::LifeGrid; sparse=true, chunklength=128, parallel=size(lg, 1)>1024)

Update `lg` one generation according to the rules of Conway's Game of Life and return it.

A Dirichlet boundary condition is applied, fixing all cells outside of the grid at zero.

`step!` runs using all available threads by default.

If `sparse` is `true`, each `$(CELLS_PER_CLUSTER)×chunklength` section of the grid is
checked for live cells so that an update of that section can be skipped if there are none.
The performance impact of this check is small but measurable, so if your grid doesn't
contain large areas devoid of living cells you should set `sparse=false`.

`chunklength` determines the size of a chunk of data that `step!` works on before proceeding
to the next chunk. 128 is chosen as the default since it strikes a good balance: it leads to
chunks large enough (5 KiB) that work isn't interrupted too often, and small enough to fit
in the L1 cache of most machines. The height of `lg` must be at least
`chunksize*Threads.nthreads()` for all threads to be fully engaged.

`parallel` determines whether `step!` will run with multiple threads. It is `true` by
default if `lg`'s height exceeds 1024, and `false` otherwise. This is a reasonable default
on most machines, but it's worth experimenting with.

# Extended Help

TODO: talk about how to specialize updatedcluster
"""
function step!(lg::LifeGrid; sparse=true, chunklength=128, parallel=size(lg, 1)>1024)
    runstep!() = if sparse
        stepraw!(lg, chunklength, Sparse())
    else
        stepraw!(lg, chunklength, :dense)
    end

    if parallel
        return runstep!()
    else
        disable_polyester_threads() do
            return runstep!()
        end
    end
end



# Type to use to tell stepraw! whether to check for dead zones while iterating
struct Sparse end



"""
    stepraw!(lg::LifeGrid, chunklength, ::Density) where Density

The back end of [`step!`](@ref); updates and returns `lg`.

If `Density` is a `LifeGame.Sparse`, the sparse algorithm is used.

# Extended help

`lg` is updated one column of clusters at a time--that is, a column of 64-bit unsigned
integers, each representing a 62-cell portion of a row (see [`updatedcluster`](@ref)).

There are 3 arrays involved in the update:

1. The grid itself, a `Matrix` of `UInt64` clusters.
2. The left buffer, a `Vector` of `UInt64`s with as many elements as the grid has rows.
3. The middle buffer, identical in type and size to the left buffer.

These are carefully updated in an order that eliminates inter-iteration dependencies when
updating a column. This allows the compiler to emit vectorized instructions and makes it
safe to operate on a single column with multiple threads. The use of only two buffer rows,
rather than an entire intermediate grid, keeps the memory footprint small and makes better
use of the cache.

Here is how one cluster on the interior of the grid is updated, following the actual order of
operations used in this function:

1. The cluster is updated to according to the life rules using the corresponding cells
   above, at the same place as, and below the cluster in the left buffer using
   `updatedcluster`:\n
   `grid[i,j] = updatedcluster(leftbuffer[i-1], leftbuffer[i], leftbuffer[i+1])`.
1. The corresponding cluster in the middle buffer has its halos updated; the corresponding
   cluster in the left buffer is used as the left value, the cluster to the right of the
   grid cluster in question as the middle value, and the cluster two to the right of the
   grid cluster in question as the right value:\n
   `rightbuffer[i] = updatedhalos(leftbuffer[i], grid[i,j+1], grid[i,j+2])`.
1. The buffers are switched in preparation for the next iteration.

Notice that no cell that is written to on this iteration is also read from; this allows for
vectorization and multithreading. Cache efficiency is increased by doing these updates in
chunks: taking a portion of a column, updating it fully, and only then moving on to the next
chunk.

Columns at the boundaries are special cases, but the order of operations is the same. The
computations are aided by having padding columns to the left and right of the first and last
active grids.
"""
function stepraw!(lg::LifeGrid{R}, chunklength, ::Density) where {R, Density}
    dense = !(Density <: Sparse)

    # Column iteration range
    Ibegin, Iend = firstindex(lg.grid, 1)+1, lastindex( lg.grid, 1)-1
    innerrange(I) = I:min(I+chunklength-1, Iend)

    # Convenience names; also helps @batch avoid allocations
    lbuf   = lg.leftcolbuffer
    mbuf   = lg.middlecolbuffer
    ldzbuf = true # only used as a temporary for a single iteration
    mdzbuf = lg.middledeadzonebuffer
    rdzbuf = lg.rightdeadzonebuffer

    # First iteration: update the halos of the second column
    firstcol  = @view lg.grid[:,begin]
    secondcol = @view lg.grid[:,begin+1]
    thirdcol  = @view lg.grid[:,begin+2]
    @batch for I in Ibegin:chunklength:Iend
        # Update the halos of the first row in preparation for inner iterations
        for i in innerrange(I)
            lbuf[i] = updatedhalos(firstcol[i], secondcol[i], thirdcol[i])
        end

        # Initialize the dead zone buffers if appropriate
        if !dense
            mdzbuf[I] = sum(@view secondcol[innerrange(I)]) != 0
            rdzbuf[I] = sum(@view thirdcol[ innerrange(I)]) != 0
        end
    end

    # Interior iterations
    @inbounds for j in firstindex(lg.grid, 2)+2:lastindex( lg.grid, 2)-1
        # Views of the current and neighboring columns
        left   = @view lg.grid[:,j-1]
        middle = @view lg.grid[:,j]
        right  = @view lg.grid[:,j+1]

        # Outer loop over chunks of rows
        @batch for I in Ibegin:chunklength:Iend
            # Swap dead zone buffers
            #if !dense
            #    ldzbuf = mdzbuf[I]
            #    mdzbuf[I] = rdzbuf[I]
            #    rdzbuf[I] = sum(@view right[innerrange(I)]) != 0
            #end

            # Update cells if appropriate
            lastI = last(innerrange(I))+1
            #if dense || ldzbuf || lbuf[I-1] != 0 || lbuf[lastI] != 0
                for i in innerrange(I)
                    left[i] = updatedcluster(lbuf[i-1], lbuf[i], lbuf[i+1], R)
                end
            #end

            # Update halos if appropriate
            #if dense || ldzbuf || mdzbuf[I] || rdzbuf[I]
                for i in innerrange(I)
                    mbuf[i] = updatedhalos(lbuf[i], middle[i], right[i])
                end
            #end
        end

        # Swap column buffers
        lbuf, mbuf = mbuf, lbuf
    end

    # Last iteration: trailing cells are zeroed
    lastcol = @view lg.grid[:,end-1]
    shift = CELLS_PER_CLUSTER - size(lg, 2)%CELLS_PER_CLUSTER + 1 # add 1 for halo
    @batch for I in Ibegin:chunklength:Iend
        # Update active cells and zero trailing cells
        if dense || ldzbuf || mdzbuf[I] || rdzbuf[I]
            for i in innerrange(I)
                updated = (updatedcluster(lbuf[i-1], lbuf[i], lbuf[i+1], R) >> shift) << shift
                lastcol[i] = updated
            end
        end
    end

    # Return the grid
    return lg
end