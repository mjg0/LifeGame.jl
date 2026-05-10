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
Rule{0x03}()
```
"""
struct Rule{NeighborSums}
    function Rule(n...)
        onbits = zero(UInt8)
        for i in n
            if i in 1:8
                onbits |= 0x01 << (i - 1)
            else
                throw(ArgumentError("Invalid rule neighbor sum $i; sums must be in [1, 8]"))
            end
        end
        return new{onbits}()
    end
end



"""
    LifeRule{Birth, Survival}

A struct holding birth and survival [`Rule`](@ref)s.

---

    LifeRule(b, s)

Return a `LifeRule` with birth rule `Rule(b...)` and survival rule `Rule(s...)`.
"""
struct LifeRule{Birth,Survival}
    function LifeRule(rule::AbstractString)
        rulematch = match(r"^B(\d*)/S(\d*)$", rule)
        if rulematch === nothing
            throw(ArgumentError("Invalid rule '$rule' supplied"))
        end
        birthnumbers, survivalnumbers = rulematch.captures
        return new{
            Rule((parse(Int, c) for c in birthnumbers)...),
            Rule((parse(Int, c) for c in survivalnumbers)...),
        }()
    end

    function LifeRule(b, s)
        return new{Rule(b...),Rule(s...)}()
    end
end

function Base.show(io::IO, ::LifeRule{B,S}) where {B,S}
    rulestr(rule) = prod(["$i" for i in rulesums(rule)]; init="")
    print(io, "B$(rulestr(B))/S$(rulestr(S))")
end



# Types to store a LifeGrid's halo matrices and grid chunk buffers
"""
    Halos{H}

A type to store the halos between the columns of a `LifeGrid`'s `grid`.

There are 4 matrices, each of the same size as `grid`: right and left halos, for both the
current and next iterations. See [`step!`](@ref)
"""
mutable struct Halos{H}
    currentleft::Matrix{H}
    nextleft::Matrix{H}
    currentright::Matrix{H}
    nextright::Matrix{H}

    Halos(H, m, n) = new{H}(ntuple(_ -> zeros(H, m, n), 4)...)
end

mutable struct Buffer{C}
    # Inner vectors store a chunk's worth of cells; there are Threads.nthreads() of each
    current::Vector{C}
    next::Vector{C}

    Buffer(::Type{C}, ::Type{H}) where {C,H} = new{C}(ntuple(_ -> zeros(C, nbits(H) + 2), 2)...)
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

[`LifePattern`](@ref)s can be inserted into a `LifeGrid` via [`insert!`](@ref).

See the extended help for [`LifeGame`](@ref) for implementation details.

---

    LifeGrid(m, n; rule="B3/S23")

Return an m×n `LifeGrid` with no living cells and rule `rule`.

---

    LifeGrid(grid::AbstractMatrix; rule="B3/S23")

Return a LifeGrid with cell values defined by `grid` with rule `rule`.

True or non-zero values indicate living cells; false or zero values indicate dead cells.
"""
struct LifeGrid{LifeRule,CType,HType,Tall,Wide} <: AbstractMatrix{Bool}
    height::Int64
    width::Int64
    grid::Matrix{CType}
    halos::Halos{HType}
    buffers::Vector{Buffer{CType}}

    # The backing array and vectors are padded, with zero cells surrounding each edge
    function LifeGrid(
        m::Integer,
        n::Integer;
        rule::AbstractString="B3/S23",
        CType::Type{<:Unsigned}=smallestuint(n + 2),
        HType::Type{<:Unsigned}=smallestuint(m),
    )
        # Halos
        haloheight = cld(m, nbits(HType))
        halowidth = cld(n, nbits(CType) - 2)
        halos = Halos(HType, haloheight, halowidth)

        # Buffers
        nbuffers = min(Threads.nthreads(), halowidth)
        buffers = [Buffer(CType, HType) for _ = 1:nbuffers]

        # Clusters
        gridheight = nbits(HType) * haloheight # store extra rows for even chunk sizes
        gridwidth = halowidth
        grid = zeros(CType, gridheight, gridwidth)

        # Size parameters
        Tall = haloheight > 1
        Wide = halowidth > 1

        return new{LifeRule(rule),CType,HType,Tall,Wide}(m, n, grid, halos, buffers)
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
    #J, shift = divrem(j-1, nbits(C)-2).+1
    J = (j - 1) ÷ (nbits(C) - 2) + 1
    shift = (j - 1) % (nbits(C) - 2) + 1

    return I, J, shift
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

    # Which halo needs to be updated?
    hidx, hmask = indexhalos(lg, i, j)
    halos = shift == 1 ? lg.halos.currentleft : lg.halos.currentright

    # Update the halo if necessary
    op = ifelse(val != zero(val), x -> x | hmask, x -> x & ~hmask)
    if shift == 1
        halos[hidx] = op(halos[hidx])
    elseif shift == nbits(C) - 2
        halos[hidx] = op(halos[hidx])
    end

    return val
end



"""
    rule(lg::LifeGrid)

Return the simulation rule governing `lg`'s evolution.
"""
rule(::LifeGrid{R}) where {R} = R



"""
    rulesums(N::Integer)
    rulesums(::Rule{N})
    rulesums(::LifeRule{B, S})

Return a vector of numbers for which the given rule specifies birth or survival.

`N` is a `UInt8`, each bit corresponding to a number 1 to 8; the numbers corresponding to on
bits are returned.

If a `LifeRule` is provided, the birth and survival sums, respectively, are returned as a
tuple.
"""
rulesums(N::Integer) = [i for i = 1:8 if N >> (i - 1) & 0x01 == 0x01]

rulesums(::Rule{N}) where {N} = rulesums(N)

rulesums(::LifeRule{B,S}) where {B,S} = rulesums(B), rulesums(S)
