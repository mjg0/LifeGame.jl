"""
    LifeGame

A module with infrastructure to simulate Conway's Game of Life.

This threaded life game implementation, which was inspired by [exrok's Rust
implementation](https://github.com/exrok/game_of_life/), uses finite grids with all cells
beyond the edges of the grid fixed at zero. It's optimized for large, dense grids.

The public interface to `LifeGame` includes the types [`LifeGame`](@ref) and
[`LifePattern`](@ref), the module [`LifePatterns`](@ref), and the functions [`step!`](@ref)
and [`insert!`](@ref).

# Examples

```jldoctest
julia> using LifeGame

julia> lg = LifeGrid(4, 5);

julia> insert!(lg, 1, 1, LifePatterns.glider)
4×5 LifeGrid:
 0  1  0  0  0
 0  0  1  0  0
 1  1  1  0  0
 0  0  0  0  0

julia> step!(lg)
4×5 LifeGrid:
 0  0  0  0  0
 1  0  1  0  0
 0  1  1  0  0
 0  1  0  0  0
```

See the extended help for an overview of `LifeGame`'s implementation.

# Extended help

## `LifeGrid` implementation

The fundamental unit of a `LifeGrid` is the **cluster**, a row of 62 cells represented by a
single `UInt64`; the two extra bits at the beginning and end of the cluster are **halo
cells**, which hold the first and last cells of the clusters to the right and left,
respectively, for numerical convenience. Zero padding is applied to each edge of the backing
array, also for numerical convenience. The storage backing a 200×300 `LifeGrid` is thus a
202×7 `Matrix{UInt64}`. This means that large grids can be stored efficiently: a
100,000×100,000 cell `LifeGrid` occupies only 1.2 GiB of memory.

## `step!` implementation

Rather than updating each bit individually, `step!` uses bitwise operations (half and full 
adders, bitshifts, **AND**s and **OR**s) on entire clusters. The function used to update a
single cluster--62 cells--compiles to 37 instructions on an x86_64 test machine. See the
extended help for [`LifeGame.updatedcluster`](@ref) for details on how the cluster update
works.

## Performance

`step!` is written to compile to highly vectorized instructions, uses CPU caches
efficiently, and is parallelized. On a laptop with an AMD 7640U, it typically takes about
250 μs to `step!` a dense 10,000×10,000 `LifeGrid`, and 50 ms to `step!` a dense
100,000×100,000 `LifeGrid`. Since keeping the CPU fed is a major bottleneck when update
operations are so fast, `step!` operates faster per cell when the grid it's working on fits
in the CPU cache.
"""
module LifeGame

using Polyester

export LifeGrid, step!



const CLUSTER_TYPE = UInt64

const CELLS_PER_CLUSTER = 8*sizeof(CLUSTER_TYPE)-2 # subtract 2 for halos at cell edges

const FIRST_BIT = CLUSTER_TYPE(2) << CELLS_PER_CLUSTER

const LAST_BIT = one(CLUSTER_TYPE)



"""
    LifeGrid <: AbstractMatrix{Bool}

The type representing a Conway's Game of Life grid.

A `LifeGrid` can be updated with the [`step!`](@ref) function.

[`LifePattern`](@ref)s can be inserted into a `LifeGrid` via [`insert!`](@ref).

See the extended help for [`LifeGame`](@ref) for implementation details.

---

    LifeGrid(m, n)

Return an m×n `LifeGrid` with no living cells.

---

    LifeGrid(grid::BitMatrix)
    LifeGrid(grid::AbstractMatrix{Bool})
    LifeGrid(grid::AbstractMatrix{Number})

Return a LifeGrid with cell values defined by `grid`.

True and non-zero values indicate living cells; false and zero values indicate dead cells.
"""
mutable struct LifeGrid <: AbstractMatrix{Bool}
    width::Int64
    grid::Matrix{CLUSTER_TYPE}
    leftcolbuffer::Vector{CLUSTER_TYPE}   # used in step!
    middlecolbuffer::Vector{CLUSTER_TYPE} # used in step!
    middledeadzonebuffer::Vector{Bool}   # used in step! when sparse=true
    rightdeadzonebuffer::Vector{Bool}    # used in step! when sparse=true

    # The backing array and vectors are padded, with zero cells surrounding each edge
    function LifeGrid(m::Integer, n::Integer)
        paddedheight = m+2
        paddedpackedwidth = cld(n, CELLS_PER_CLUSTER)+2
        grid = zeros(CLUSTER_TYPE, paddedheight, paddedpackedwidth)
        buffer = zeros(CLUSTER_TYPE, paddedheight)
        return new(n, grid, buffer, deepcopy(buffer), zeros(Bool, m+2), zeros(Bool, m+2))
    end

    function LifeGrid(grid::BitArray)
        lg = LifeGrid(size(grid)...)
        lg .= grid
        return lg
    end

    LifeGrid(grid::AbstractMatrix{Bool}) = LifeGrid(BitArray(grid))

    LifeGrid(grid::AbstractMatrix{<:Number}) = LifeGrid(grid.!=0)
end



# Implement AbstractArray interface for LifeGrid
Base.size(lg::LifeGrid) = size(lg.grid, 1)-2, lg.width

function indexlifegrid(i, j)
    I = i+1 # skip the padding cells
    J = (j-1)÷CELLS_PER_CLUSTER+2 # add one for 1-index arrays
    shift = (j-1)%CELLS_PER_CLUSTER+1
    return I, J, shift
end

Base.@propagate_inbounds function Base.getindex(lg::LifeGrid, i::Integer, j::Integer)
    I, J, shift = indexlifegrid(i, j)
    return ((lg.grid[I,J] << shift) & FIRST_BIT) == FIRST_BIT 
end

Base.@propagate_inbounds function Base.setindex!(lg::LifeGrid, val::Bool,
                                                 i::Integer, j::Integer)
    I, J, shift = indexlifegrid(i, j)
    cluster = lg.grid[I,J]
    lg.grid[I,J] = ifelse(val,
                          cluster |   FIRST_BIT >> shift,
                          cluster & ~(FIRST_BIT >> shift))
    return val
end

Base.@propagate_inbounds function Base.setindex!(lg::LifeGrid, val::Number,
                                                 i::Integer, j::Integer)
    return lg[i,j] = val != zero(typeof(val))
end



"""
    updatedcluster(above::Integer, current::Integer, below::Integer)

Update and return `current` given adjacent clusters `above` and `below` using bit-twiddling.

`updatedcluster` uses bitwise operations to efficiently count the number of neighbors of
each cell in the `current` cluster and return the updated cluster. It employs half and full
adders to calculate the sums and carries needed.

# Extended help

It's helpful to think of operations on clusters (`UInt64`s representing 62-cell portions of
rows) happening to a single bit, but since all the operators used are bitwise, any
computation occurs simultaneously on the entire cluster. Some of the explanations that
follow will reference a single bit (or cell) for conceptual simplicity.

At a high level, `updatedcluster` sums the neighboring cells of the cell in question using
half and full adders. Since these are bitwise operators, each bit in the sum must be
represented separately: the ones bit, the twos bit, and the fours bit of the binary number
are all calculated. Once these sums are obtained, they are used to determine whether the
cell will be alive in the next generation.

The algorithm used by `updatedcluster` is as follows:

1. Determine the sum of the middle column: a full adder is used to sum `above`, `middle`, \
   and `below`. This results in a "sum" (the bit representing the ones place) and a \
   "remainder" (the bit representing the twos place).
2. Determine the sums of the left and right columns: these are obtained by shifting the \
   bits of the middle one to the left and one to the right, repsectively.
3. Determine the sum of the cells directly above and below the cell in question: a half \
   adder is used to sum `above` and `below`.
4. Sum the ones place from the left and right columns and above/below: a full adder is \
   used to get the ones and twos place from these three previous sums.
5. Sum the twos place from the left and right columns and above/below: a full adder is \
   used to get the twos and threes place from these three previous remainders.
6. Sum the twos places: a half adder is used to add the remainder from **4** and the sum \
   from **5**.
7. Update the cluster: `current` is modified three times before being returned:
    - All cells for which the first bit from **4** are flipped on if they aren't on \
      already. This corresponds to odd sums; odd sums that aren't 3 are turned off in the
      next steps.
    - All cells for which the second bit is off are turned off. No cell for which the twos \
      bit sum isn't on can be alive, since both 2 (0x10) and 3 (0x11) have the twos bit on.
    - All cells for which the partial third bit sum from **5** is on are turned off. If \
      this bit is on, the sum of the cell's neighbors is at least 4.
"""
function updatedcluster(above, current, below)
    # Bitwise half and full adders, returning a sum and a remainder
    halfadder(x, y)    = x ⊻ y,     x & y
    fulladder(x, y, z) = x ⊻ y ⊻ z, x&y | x&z | y&z

    # Sums and remainders of each column
    middlesum, middlerem = halfadder(above, below) # excludes the middle cell
    basesum,   baserem   = fulladder(above, current, below)
    leftsum,   leftrem   = basesum << 1, baserem << 1
    rightsum,  rightrem  = basesum >> 1, baserem >> 1

    # Sums in each bit: the first bit representing 1, the second 2, the third 4, etc.
    bit1,  bit2a = fulladder(leftsum, middlesum, rightsum)
    bit2b, bit3a = fulladder(leftrem, middlerem, rightrem)
    bit2,  bit3b = halfadder(bit2a, bit2b)
    # bit3b is unnecessary since if bit3a is on, that already implies too many neighbors;
    # there is thus also no need to determine the value of the fourth bit

    # Update and return current
    current |=  bit1  # Switch on cells with an odd neighbor count (i.e. the one bit on)
    current &=  bit2  # ...but only keep cells on if the two bit is also on
    current &= ~bit3a # ...and if there are less than 4 living neighbors
    return current
end



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
function stepraw!(lg::LifeGrid, chunklength, ::Density) where Density
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
            if !dense
                ldzbuf = mdzbuf[I]
                mdzbuf[I] = rdzbuf[I]
                rdzbuf[I] = sum(@view right[innerrange(I)]) != 0
            end

            # Update cells if appropriate
            lastI = min(I+chunklength-1, lastindex(lbuf))
            if dense || ldzbuf || lbuf[I-1] != 0 || lbuf[lastI] != 0
                for i in innerrange(I)
                    left[i] = updatedcluster(lbuf[i-1], lbuf[i], lbuf[i+1])
                end
            end

            # Update halos if appropriate
            if dense || ldzbuf || mdzbuf[I] || rdzbuf[I]
                for i in innerrange(I)
                    mbuf[i] = updatedhalos(lbuf[i], middle[i], right[i])
                end
            end
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
                updated = (updatedcluster(lbuf[i-1], lbuf[i], lbuf[i+1]) >> shift) << shift
                lastcol[i] = updated
            end
        end
    end

    # Return the grid
    return lg
end



include("LifePattern.jl")

include("LifePatterns.jl")



end # module