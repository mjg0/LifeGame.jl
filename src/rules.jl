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
