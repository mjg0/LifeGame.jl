export LifeGrid, rule



"""
    Rule{NeighborSums}

A birth or survival rule, where `NeighborSums` stores the numbers for the rule.

`NeighborSums` is an 8-bit unsigned integer, where the `n`th bit being on means that number
leads to birth or survival. For example, high life, with the rule B36/S23, would have birth
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
            if i in n
                onbits |= 0x01<<(i-1)
            else
                throw(ArgumentError("Invalid rule neighbor sum $i; sums must be in [1, 8]"))
            end
        end
        return new{onbits}()
    end
end



"""
    LifeRule{Birth,Survival}

A struct holding birth and survival [`Rule`](@ref)s.

---

    LifeRule(b, s)

Return a `LifeRule` with birth rule `Rule(b...)` and survival rule `Rule(s...)`.
"""
struct LifeRule{Birth, Survival}
    function LifeRule(rule::AbstractString)
        rulematch = match(r"^B(\d*)/S(\d*)$", rule)
        if rulematch === nothing
            throw(ArgumentError("Invalid rule '$rule' supplied"))
        end
        birthnumbers, survivalnumbers = rulematch.captures
        return new{Rule((parse(Int, c) for c in birthnumbers   )...),
                   Rule((parse(Int, c) for c in survivalnumbers)...)}()
    end

    function LifeRule(b, s)
        return new{Rule(b...), Rule(s...)}()
    end
end

function Base.show(io::IO, ::LifeRule{B, S}) where {B, S}
    rulestr(rule) = prod(["$i" for i in rulesums(rule)]; init="")
    print(io, "B$(rulestr(B))/S$(rulestr(S))")
end



"""
    LifeGrid <: AbstractMatrix{Bool}

The type representing a grid for the simulation of 2-D cellular automata.

Each `LifeGrid` has an associated rule for the automata that determines how many living
neighbor cells lead to birth and how many to survival. It's formatted as `Bm.../Sn...`,
where `m...` are the non-deliminted numbers leading to birth, and `n...` are the
non-deliminted numbers leading to survival. The default is the rule for Conway's Game of
Life, `B3/S23`. Query the rule for a given `LifeGrid` with the [`rule`](@ref) function.

A `LifeGrid` can be updated with the [`step!`](@ref) function.

[`LifePattern`](@ref)s can be inserted into a `LifeGrid` via [`insert!`](@ref).

See the extended help for [`LifeGame`](@ref) for implementation details.

---

    LifeGrid(m, n; rule="B3/S23")

Return an m×n `LifeGrid` with no living cells and rule `rule`.

---

    LifeGrid(grid::BitMatrix; rule="B3/S23")
    LifeGrid(grid::AbstractMatrix{Bool}; rule="B3/S23")
    LifeGrid(grid::AbstractMatrix{Number}; rule="B3/S23")

Return a LifeGrid with cell values defined by `grid` with rule `rule`.

True and non-zero values indicate living cells; false and zero values indicate dead cells.
"""
mutable struct LifeGrid{LifeRule} <: AbstractMatrix{Bool}
    width::Int64
    grid::Matrix{CLUSTER_TYPE}
    leftcolbuffer::Vector{CLUSTER_TYPE}   # used in step!
    middlecolbuffer::Vector{CLUSTER_TYPE} # used in step!

    # The backing array and vectors are padded, with zero cells surrounding each edge
    function LifeGrid(m::Integer, n::Integer; rule::AbstractString="B3/S23")
        # Create grid based on the supplied size
        paddedheight = m+2
        paddedpackedwidth = cld(n, CELLS_PER_CLUSTER)+2
        grid = zeros(CLUSTER_TYPE, paddedheight, paddedpackedwidth)

        # Buffers
        lcbuf = zeros(CLUSTER_TYPE, paddedheight)
        mcbuf = deepcopy(lcbuf)

        # Return the LifeGrid
        return new{LifeRule(rule)}(n, grid, lcbuf, mcbuf)
    end

    function LifeGrid(grid::BitArray; kw...)
        lg = LifeGrid(size(grid)...; kw...)
        lg .= grid
        return lg
    end

    LifeGrid(grid::AbstractMatrix{Bool}; kw...) = LifeGrid(BitArray(grid); kw...)

    LifeGrid(grid::AbstractMatrix{<:Number}; kw...) = LifeGrid(grid.!=0; kw...)
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
    rule(lg::LifeGrid)

Return the simulation rule governing `lg`'s evolution.
"""
rule(::LifeGrid{R}) where R = R


rulesums(N) = [i for i in 1:8 if N>>(i-1)&0x01==0x01]
rulesums(::Rule{N}) where N = rulesums(N)
rulesums(::LifeRule{B, S}) where {B, S} = rulesums(B), rulesums(S)
