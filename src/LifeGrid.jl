export LifeGrid, rule



"""
    LifeGrid <: AbstractMatrix{Bool}

The type representing a grid for the simulation of 2-D cellular automata.

Each `LifeGrid` has an associated rule for the automata that determines how many living
neighbor cells lead to birth and how many to survival. It's formatted as `Bm.../Sn...`,
where `m...` are the non-deliminted numbers leading to birth, and `n...` are the
non-deliminted numbers leading to survival. The default is the rule for Conway's Game of
Life, `B3/S23`. Query the rule for a given `LifeGrid` with the [`rule`](@ref) function.

A `LifeGrid` can be advanced one generation with the [`step!`](@ref) function.

[`LifePattern`](@ref)s and other matrices can be inserted into a `LifeGrid` with
[`insert!`](@ref).

See the extended help for [`LifeGame`](@ref) for implementation details.

---

    LifeGrid(m, n; rule="B3/S23", CType, HType)

Return an m×n `LifeGrid` with no living cells and rule `rule`.

`CType` determines the type used to store clusters of cells, which are N-2 cells packed into
unsigned integers with N bits. It defaults to the smallest unsigned type that can fit the
whole grid's width, `UInt64` if the grid is 31 or more cells wide. Since entire clusters are
updated with bitwise operations even if some of the trailing bits are empty, the grid widths
that make [`step!`](@ref) most efficient are `8*2ⁿ-2, n∈[0,3]`: 6, 14, 30, and multiples of
62.

Likewise, `HType` determines the type in which halos are stored, with similar defaults, up
to `UInt64` if the grid is more than 32 bits tall. Since entire "chunks" (sections of each
column with as many elements as `HType` has bits) are updated even if some trailing clusters
are empty, grid heights of 8, 16, 32, and multiples of 64 will get the most performance out
of [`step!`](@ref).

---

    LifeGrid(grid::AbstractMatrix; rule="B3/S23", CType, HType)

Return a LifeGrid with cell values defined by `grid` with rule `rule`.

True or non-zero values indicate living cells; false or zero values indicate dead cells.
"""
struct LifeGrid{R,CType,HType,Tall,Wide} <: AbstractMatrix{Bool}
    # See the help message for `updatecolumn!` for an explanation of these members.
    height::Int64
    width::Int64
    grid::Matrix{CType}
    halos::Halos{HType}
    buffers::Vector{Vector{CType}}
    changed::SparseActiveFlags

    function LifeGrid(
        m::Integer,
        n::Integer;
        rule = "B3/S23",
        CType::Type{<:Unsigned} = smallestuint(n + 2),
        HType::Type{<:Unsigned} = smallestuint(m),
    )
        # Normalize rule
        R = rule isa LifeRule ? rule : LifeRule(rule)

        # Halos
        haloheight = cld(m, nbits(HType))
        halowidth = cld(n, nbits(CType) - 2)
        halos = Halos(HType, haloheight, halowidth)

        # Buffers
        nbuffers = 2 * min(Threads.nthreads(), halowidth)
        buffers = [zeros(CType, nbits(HType) + 2) for _ = 1:nbuffers]

        # Clusters
        gridheight = nbits(HType) * haloheight # store extra rows for even chunk sizes
        gridwidth = halowidth
        grid = zeros(CType, gridheight, gridwidth)

        # Sparse update state
        changed = SparseActiveFlags(haloheight, halowidth)

        # Size parameters
        Tall = haloheight > 1
        Wide = halowidth > 1

        return new{typeof(R),CType,HType,Tall,Wide}(
            m,
            n,
            grid,
            halos,
            buffers,
            changed,
        )
    end

    function LifeGrid(grid::AbstractMatrix{T}; kw...) where {T<:Number}
        lg = LifeGrid(size(grid)...; kw...)
        lg .= grid .!= zero(T)
        return lg
    end
end



"""
    rule(lg::LifeGrid)

Return the simulation rule governing `lg`'s evolution.
"""
rule(::LifeGrid{R,C,H,Tall,Wide}) where {R,C,H,Tall,Wide} = R()
