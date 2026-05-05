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
    LifeRule{Birth, Survival}

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

A `LifeGrid` can be advanced one generation with the [`step!`](@ref) function.

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
mutable struct LifeGrid{LifeRule,CType,HType} <: AbstractMatrix{Bool}
    # TODO: parameterize, LifeGrid{T,H,LifeRule}
    height::Int64
    width::Int64
    grid::Matrix{CType}
    colbuffers1::Vector{Vector{CType}}
    colbuffers2::Vector{Vector{CType}}
    inhalosleft::Matrix{HType}
    inhalosright::Matrix{HType}
    outhalosleft::Matrix{HType}
    outhalosright::Matrix{HType}

    # The backing array and vectors are padded, with zero cells surrounding each edge
    function LifeGrid(m::Integer, n::Integer; rule::AbstractString="B3/S23")
        uint(N) = if N<=8
            UInt8
        elseif N<=16
            UInt16
        elseif N<=32
            UInt32
        else
            UInt64
        end

        HType = UInt64 # TODO uint(m)
        CType = UInt64 # TODO uint(n+2)
        cellspercluster = 8*sizeof(CType)-2

        # Buffers
        bufferheight = cld(m+2, 8*sizeof(HType))
        bufferwidth = cld(n, cellspercluster)+2
        halos = ntuple(_->zeros(HType, bufferheight, bufferwidth), 4)
        colbuffers = ntuple(_->[zeros(CType, 8*sizeof(HType)+2) for _ in 1:Threads.nthreads()], 2)

        # Grid
        gridheight = 8*sizeof(HType)*bufferheight+2
        gridwidth = bufferwidth
        grid = zeros(CType, gridheight, gridwidth)

        # Return the LifeGrid
        return new{LifeRule(rule),CType,HType}(m, n, grid, colbuffers..., halos...)
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
Base.size(lg::LifeGrid) = lg.height, lg.width

function indexlifegrid(i, j)
    I = i+1 # skip the padding cells
    J = (j-1)÷CELLS_PER_CLUSTER+2 # add one for 1-index arrays
    shift = (j-1)%CELLS_PER_CLUSTER+1
    return I, J, shift
end

Base.@propagate_inbounds function Base.getindex(lg::LifeGrid{R,C,H}, i::Integer, j::Integer) where {R,C,H}
    firstbit = one(C)<<(8*sizeof(C)-1)
    I, J, shift = indexlifegrid(i, j)
    return ((lg.grid[I,J] << shift) & firstbit) == firstbit
end

Base.@propagate_inbounds function
Base.setindex!(lg::LifeGrid{R,C,H}, val::Number, i::Integer, j::Integer) where {R,C,H}
    Cbits = 8*sizeof(C)
    Hbits = 8*sizeof(H)

    firstbit = one(C) << (Cbits - 1)
    I, J, shift = indexlifegrid(i, j)

    cellmask = firstbit >> shift

    cluster = lg.grid[I, J]
    lg.grid[I, J] = ifelse(val != zero(val),
                           cluster | cellmask,
                           cluster & ~cellmask)

    # Which vertical halo chunk contains row i?
    hI = (I - 2) ÷ Hbits + 1
    k  = (I - 2) % Hbits + 1
    hmask = one(H) << (Hbits - k)

    # First active bit in a cluster: left edge.
    if shift == 1
        halo = lg.inhalosleft[hI, J]
        lg.inhalosleft[hI, J] = ifelse(val != zero(val),
                                       halo | hmask,
                                       halo & ~hmask)
    end

    # Last active bit in a cluster: right edge.
    if shift == Cbits-2
        halo = lg.inhalosright[hI, J]
        lg.inhalosright[hI, J] = ifelse(val != zero(val),
                                        halo | hmask,
                                        halo & ~hmask)
    end

    return val
end



"""
    rule(lg::LifeGrid)

Return the simulation rule governing `lg`'s evolution.
"""
rule(::LifeGrid{R}) where R = R



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
rulesums(N::Integer) = [i for i in 1:8 if N>>(i-1)&0x01==0x01]

rulesums(::Rule{N}) where N = rulesums(N)

rulesums(::LifeRule{B, S}) where {B, S} = rulesums(B), rulesums(S)
