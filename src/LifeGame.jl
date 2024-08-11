"""
    LifeGame

A module for simulating Conway's Game of Life and similar cellular automata.

This implementation uses finite grids with all cells beyond the edges of the grid fixed at
zero. It's optimized for large, dense grids and has impressive performance.

The public interface to `LifeGame` includes the types [`LifeGame`](@ref) and
[`LifePattern`](@ref), the module [`LifePatterns`](@ref), and the functions [`step!`](@ref),
[`insert!`](@ref), and [`rule`](@ref).

# Examples

```jldoctest
julia> using LifeGame

julia> lg = LifeGrid(4, 5, rule="B36/S23"); # simulate highlife

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
single cluster--62 cells--compiles to 40 instructions for Conway's Game of Life on an x86_64
test machine. See the extended help for [`LifeGame.updatedcluster`](@ref) for details on how
the cluster update works and how to specialize a rule to improve performance.

## Performance

`step!` is written to compile to highly vectorized instructions, uses CPU caches
efficiently, and is parallelized. On a laptop with an AMD 7640U, it typically takes about
300 μs to `step!` a 10,000×10,000 `LifeGrid` and 50 ms to `step!` a 100,000×100,000
`LifeGrid` using the Conway's Game of Life rule. Since keeping the CPU fed is a major
bottleneck when update operations are so fast, `step!` operates faster per cell when the
grid it's working on fits in the CPU cache.
"""
module LifeGame



using Polyester



const CLUSTER_TYPE = UInt64

const CELLS_PER_CLUSTER = 8*sizeof(CLUSTER_TYPE)-2 # subtract 2 for halos at cell edges

const FIRST_BIT = CLUSTER_TYPE(2)<<CELLS_PER_CLUSTER # 0b100...000

const LAST_BIT = one(CLUSTER_TYPE)                   # 0b000...001

const DEFAULT_CHUNK_SIZE = 64



include("LifeGrid.jl")

include("updatedcluster.jl")

include("step.jl")

include("LifePattern.jl")

include("LifePatterns.jl")



end # module