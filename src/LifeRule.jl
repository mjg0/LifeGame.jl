"""
    Rule(n::Integer...)

Return a rule mask for which neighbor sums in `n` lead to birth or survival.

The returned mask is a 16-bit unsigned integer. If the `n`th bit (counting from zero) is on,
having `n` neighbors leads to birth or survival. For example, highlife, with the rule
"B36/S23", has birth rule `0b0000000001001000` and survival rule
`0b0000000000001100`.

# Examples

```jldoctest
julia> Rule(1, 2)
0x0006
```
"""
function Rule(n...)
    onbits = zero(UInt16)
    for i in n
        if i in 0:8
            onbits |= 0x0001 << i
        else
            throw(ArgumentError("Invalid rule neighbor sum $i; sums must be in [0, 8]"))
        end
    end
    return onbits
end



"""
    LifeRule{Birth, Survival}

A struct holding birth and survival [`Rule`](@ref)s.

---

    LifeRule(b, s)

Return a `LifeRule` with birth rule mask `Rule(b...)` and survival rule mask `Rule(s...)`.
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

Return a vector of numbers for which the given rule specifies birth or survival.

`N` is a `UInt16`, each bit corresponding to a number 0 to 8; the numbers corresponding to
on bits are returned.

If a `LifeRule` is provided, the birth and survival sums, respectively, are returned as a
tuple.
"""
rulesums(N::Integer) = [i for i in 0:8 if (N >> i) & 0x01 == 0x01]

rulesums(::LifeRule{B,S}) where {B,S} = rulesums(B), rulesums(S)
