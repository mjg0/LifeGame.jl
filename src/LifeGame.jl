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
julia> lg = LifeGrid(4, 5, rule="B36/S23");

julia> insert!(lg, 1, 1, LifePatterns.glider1)
4×5 LifeGrid{B36/S23, UInt8, UInt8}:
 0  1  0  0  0
 0  0  1  0  0
 1  1  1  0  0
 0  0  0  0  0

julia> step!(lg)
4×5 LifeGrid{B36/S23, UInt8, UInt8}:
 0  0  0  0  0
 1  0  1  0  0
 0  1  1  0  0
 0  1  0  0  0
```

See the extended help for an overview of `LifeGame`'s implementation.

# Extended help

## `LifeGrid` implementation

The fundamental unit of a `LifeGrid` is the **cluster**, usually a row of 62 cells
represented by a single `UInt64`; the two extra bits at the beginning and end of the cluster
are **halo cells**, which hold the first and last cells of the clusters to the right and
left, respectively, for numerical convenience. Clusters can be smaller unsigned integers;
for example, a 6x6 `LifeGrid` has clusters of type `UInt8` by default.

Inter-iteration halos are stored in packed matrices. A *chunk* of a grid is a portion of a
column of clusters of height equal to the number of bits in the each element of these
matrices. Whole chunks are updated simultaneously, so most grids will have some extra rows
to ensure a height the number of bits in the halo type divides evenly.

As an example, the storage backing a 200×300 `LifeGrid` is thus a 256×5 `Matrix{UInt64}`,
plus a few much smaller matrices. This means that large grids can be stored efficiently:
a 100,000×100,000 cell `LifeGrid` occupies 1.3 GiB of memory.

## `step!` implementation

Rather than updating each bit individually, `step!` uses bitwise operations (half and full 
adders, bitshifts, **AND**s and **OR**s) on entire clusters. The function used to update a
single cluster--62 cells--compiles to 40 instructions for Conway's Game of Life on an x86_64
test machine. See the extended help for [`LifeGame.updatedcluster`](@ref) for details on how
the cluster update works and how to specialize a rule to improve performance.

## Performance

`step!` is written to compile to highly vectorized instructions, uses CPU caches
efficiently, and is parallelized. On a laptop with an AMD 7640U, it typically takes about
4 μs to `step!` a 1,000×1,000 `LifeGrid` and 50 ms to `step!` a 100,000×100,000 `LifeGrid`
when using the Conway's Game of Life rule. Since keeping the CPU fed is a major bottleneck
when update operations are so fast, `step!` is usually memory-bound for larger grids on
modern hardware.
"""
module LifeGame



using Polyester



include("utils.jl")

include("LifeGrid.jl")

include("BitReaderWriter.jl")

include("io.jl")

include("updatedcluster.jl")

include("step.jl")

include("LifePattern.jl")

include("LifePatterns.jl")



end # module
