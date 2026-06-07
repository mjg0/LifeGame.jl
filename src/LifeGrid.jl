export LifeGrid, rule



"""
    Rule{NeighborSums}

A birth or survival rule, where `NeighborSums` stores the numbers for the rule.

`NeighborSums` is an 8-bit unsigned integer, where the `n`th bit being on means that number
leads to birth or survival. For example, highlife, with the rule B36/S23, would have birth
rule `Rule{0b00100100}` and survival rule `Rule{0b00000110}`.

---

    Rule(n::Integer...)

Return a `Rule` for which neighbor sums in `n` lead to birth or survival.

# Examples

```jldoctest
julia> Rule(1, 2)
Rule{0x03}
```
"""
struct Rule{NeighborSums} end

function Rule(n...)
    onbits = zero(UInt8)
    for i in n
        if i in 1:8
            onbits |= 0x01 << (i - 1)
        else
            throw(ArgumentError("Invalid rule neighbor sum $i; sums must be in [1, 8]"))
        end
    end
    return Rule{onbits}
end



"""
    LifeRule{Birth, Survival}

A struct holding birth and survival [`Rule`](@ref)s.

---

    LifeRule(b, s)

Return a `LifeRule` with birth rule `Rule(b...)` and survival rule `Rule(s...)`.
"""
struct LifeRule{Birth,Survival} end

function LifeRule(rule::AbstractString)
    rulematch = match(r"^B(\d*)/S(\d*)$", rule)
    if rulematch === nothing
        throw(ArgumentError("Invalid rule '$rule' supplied"))
    end
    birthnumbers, survivalnumbers = rulematch.captures
    return LifeRule{
        Rule((parse(Int, c) for c in birthnumbers)...),
        Rule((parse(Int, c) for c in survivalnumbers)...),
    }()
end

function LifeRule(b, s)
    return LifeRule{Rule(b...),Rule(s...)}()
end

function Base.show(io::IO, ::LifeRule{B,S}) where {B,S}
    rulestr(rule) = prod(["$i" for i in rulesums(rule)]; init = "")
    print(io, "B$(rulestr(B))/S$(rulestr(S))")
end



# Types to store a LifeGrid's halo matrices and grid chunk buffers
"""
    Halos{H}

A type to store the halos between the columns of a `LifeGrid`'s `grid`.

There are 4 matrices, each of the same size as `grid`: right and left halos, for both the
current and next iterations.
"""
mutable struct Halos{H}
    currentleft::Matrix{H}
    nextleft::Matrix{H}
    currentright::Matrix{H}
    nextright::Matrix{H}

    Halos(H, m, n) = new{H}(ntuple(_ -> zeros(H, m, n), 4)...)
end



"""
    smallestuint(N)

Return the smallest `Unsigned` type that has at least `min(N, 64)` bits.
"""
function smallestuint(N)
    return if N ≤ 8
        UInt8
    elseif N ≤ 16
        UInt16
    elseif N ≤ 32
        UInt32
    else
        UInt64
    end
end



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

    function LifeGrid(
        m::Integer,
        n::Integer;
        rule = "B3/S23",
        CType::Type{<:Unsigned} = smallestuint(n + 2),
        HType::Type{<:Unsigned} = smallestuint(m),
    )
        # Normalize rule
        rule = rule isa LifeRule ? rule : LifeRule(rule)

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

        # Size parameters
        Tall = haloheight > 1
        Wide = halowidth > 1

        return new{typeof(rule),CType,HType,Tall,Wide}(m, n, grid, halos, buffers)
    end

    function LifeGrid(grid::AbstractMatrix{T}; kw...) where {T<:Number}
        lg = LifeGrid(size(grid)...; kw...)
        lg .= grid .!= zero(T)
        return lg
    end
end



"""
    indexlifegrid(lg::LifeGrid, i, j)

Return the coordinates and a shift translating index `(i, j)` to a cell in `lg.grid`.

Coordinates are two indices and the shift is a number up to `8*sizeof(eltype(lg.grid))`.

For coordinates and shift obtained thus:

```julia
I, J, shift = indexlifegrid(lg, i, j)
```

...`lg[i,j]` is true only if `(lg.grid[I,J] << shift` has its highest bit on.
"""
Base.@propagate_inbounds @inline function indexlifegrid(
    ::LifeGrid{R,C,H},
    i,
    j,
) where {R,C,H}
    I = i
    J = (j - 1) ÷ (nbits(C) - 2) + 1
    shift = (j - 1) % (nbits(C) - 2) + 1

    return I, J, shift
end



"""
    gridchunk(lg::LifeGrid, I, J)

Return a view of the `(I, J)`th chunk of `lg`.

A "chunk" is a section of a column as many elements long as `lg`'s halo type has bits. This
is usually 64 clusters for large grids, but could be as few as 8.
"""
Base.@propagate_inbounds @inline function gridchunk(
    lg::LifeGrid{R,C,H,Tall,Wide},
    I,
    J,
) where {R,C,H,Tall,Wide}
    i = (I - 1) * nbits(H) + 1
    return view(lg.grid, i:(i+nbits(H)-1), J)
end



"""
    indexhalos(lg::LifeGrid, i, j)

Return an index and mask for indexing the halos in `lg`, given conceptual index `(i, j)`.

Use it thus:

```julia
halos = lg.halos.leftnext # output left halos
I, mask = indexhalos(lg, i, j)

# Set to true
halos[I] = halos |  mask

# Set to false
halos[I] = halos & ~mask
```
"""
Base.@propagate_inbounds @inline function indexhalos(
    lg::LifeGrid{R,C,H},
    i::Integer,
    j::Integer,
) where {R,C,H}
    I, J, _ = indexlifegrid(lg, i, j)

    # The index in the halo array
    hidx = CartesianIndex((I - 1) ÷ nbits(H) + 1, J)

    # A single-bit mask with the appropriate shift
    k = (I - 1) % nbits(H) + 1
    mask = one(H) << (nbits(H) - k)

    return hidx, mask
end



"""
    syncclusterhalos!(lg::LifeGrid, I, J)

Update `lg`'s current halo matrices for the `(I, J)`th cluster of `lg.grid`.

This should be called after directly mutating `lg.grid[I, J]` without using `setindex!`. It
copies the first and last real cell bits of the cluster into `lg.halos.currentleft` and
`lg.halos.currentright`.
"""
Base.@propagate_inbounds function syncclusterhalos!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    I::Integer,
    J::Integer,
) where {R,C,H,Tall,Wide}
    hidx = CartesianIndex((I - 1) ÷ nbits(H) + 1, J)
    k = (I - 1) % nbits(H) + 1
    hmask = one(H) << (nbits(H) - k)

    cluster = lg.grid[I, J]

    leftbit = (cluster >> (nbits(C) - 2)) & one(C)
    rightbit = (cluster >> 1) & one(C)

    lg.halos.currentleft[hidx] = ifelse(
        leftbit == one(C),
        lg.halos.currentleft[hidx] | hmask,
        lg.halos.currentleft[hidx] & ~hmask,
    )

    lg.halos.currentright[hidx] = ifelse(
        rightbit == one(C),
        lg.halos.currentright[hidx] | hmask,
        lg.halos.currentright[hidx] & ~hmask,
    )

    return nothing
end



# AbstractArray interface for LifeGrid
Base.size(lg::LifeGrid) = lg.height, lg.width

Base.@propagate_inbounds function Base.getindex(
    lg::LifeGrid{R,C,H},
    i::Integer,
    j::Integer,
) where {R,C,H}
    I, J, shift = indexlifegrid(lg, i, j)

    return ((lg.grid[I, J] << shift) & highbit(C)) == highbit(C)
end

Base.@propagate_inbounds function Base.setindex!(
    lg::LifeGrid{R,C,H},
    val::Number,
    i::Integer,
    j::Integer,
) where {R,C,H}
    I, J, shift = indexlifegrid(lg, i, j)

    # Update the cluster
    cellmask = highbit(C) >> shift
    cluster = lg.grid[I, J]
    lg.grid[I, J] = ifelse(val != zero(val), cluster | cellmask, cluster & ~cellmask)

    # Update halos if needed
    if shift == 1 || shift == nbits(C) - 2
        syncclusterhalos!(lg, I, J)
    end

    return val
end



"""
    rule(lg::LifeGrid)

Return the simulation rule governing `lg`'s evolution.
"""
rule(::LifeGrid{R,C,H,Tall,Wide}) where {R,C,H,Tall,Wide} = R()



"""
    rulesums(N::Integer)
    rulesums(::Type{Rule{N}})
    rulesums(::Type{LifeRule{B, S}})

Return a vector of numbers for which the given rule specifies birth or survival.

`N` is a `UInt8`, each bit corresponding to a number 1 to 8; the numbers corresponding to on
bits are returned.

If a `LifeRule` is provided, the birth and survival sums, respectively, are returned as a
tuple.
"""
rulesums(N::Integer) = [i for i = 1:8 if N >> (i - 1) & 0x01 == 0x01]

rulesums(::Rule{N}) where {N} = rulesums(N)
rulesums(::Type{R}) where {R<:Rule} = rulesums(R())

rulesums(::LifeRule{B,S}) where {B,S} = rulesums(B), rulesums(S)
rulesums(::Type{R}) where {R<:LifeRule} = rulesums(R())
