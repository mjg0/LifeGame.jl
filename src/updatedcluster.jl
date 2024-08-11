"""
    bitsums(above, current, below)

Return a 4-tuple containing the bitwise sum of the neighbor bits of the cluster `current`.

The first element of the returned tuple represents the ones place, the second the twos
place, the third the fours place, and the fourth the eights place.

# Extended Help

As an example, given these `above`, `current`, and `below`:

```julia
above   = 0b00100000    # 0 0 1 0
current = 0b01110000    # 0 1 1 1
below   = 0b01100000    # 0 1 1 0
```

...the returned bitwise sums would be:

```julia
bit1    = 0b00010000    # 0 0 1 1   ones place
bit2    = 0b10010000    # 1 0 0 1   twos place
bit3    = 0b01100000    # 0 1 1 0   fours place
bit4    = 0b00000000    # 0 0 0 0   eights place
                        # 2 4 5 3   TOTALS
```

...meaning that the first bit of `current` has 2 neighbors, the second has 4, the third has
5, and the fourth has 3.
"""
@inline function bitsums(above, current, below)
    # Bitwise half and full adders, returning a sum and a remainder
    halfadder(x, y)    = x ⊻ y,     x & y
    fulladder(x, y, z) = x ⊻ y ⊻ z, x&y | x&z | y&z

    # Sums and remainders of each column
    middlesum, middlerem = halfadder(above, below) # excludes the middle cell
    basesum,   baserem   = fulladder(above, current, below)
    leftsum,   leftrem   = basesum << 1, baserem << 1
    rightsum,  rightrem  = basesum >> 1, baserem >> 1

    # Sums in each bit: the 1st bit represents 1, the 2nd 2, the 3rd 4, and the 4th 8
    bit1,  bit2a = fulladder(leftsum, middlesum, rightsum)
    bit2b, bit3a = fulladder(leftrem, middlerem, rightrem)
    bit2,  bit3b = halfadder(bit2a, bit2b)
    bit3,  bit4  = halfadder(bit3a, bit3b)

    # Return all 4 bits
    return bit1, bit2, bit3, bit4
end



"""
    alivebits(bit1, bit2, bit3, bit4, ::Rule)

Return a bitmask indicating which cells will be alive in the next generation.

The bits that remain alive are determined according to the bitwise sum of neighbor cells
represented by `bit*` (see [`bitsums`](@ref)) and the provided [`Rule`](@ref).
"""
@generated function alivebits(bit1::T, bit2::T, bit3::T, bit4::T, ::Rule{N}) where {T, N}
    # "On" bits for possible sums of 1 through 8
    onbits = (:( bit1 & ~bit2 & ~bit3 & ~bit4), # 1: 1st bit
              :(~bit1 &  bit2 & ~bit3 & ~bit4), # 2: 2nd bit
              :( bit1 &  bit2 & ~bit3 & ~bit4), # 3: 1st, 2nd bits
              :(~bit1 & ~bit2 &  bit3 & ~bit4), # 4: 3rd bit
              :( bit1 & ~bit2 &  bit3 & ~bit4), # 5: 1st, 3rd bits
              :(~bit1 &  bit2 &  bit3 & ~bit4), # 6: 2nd, 3rd bits
              :( bit1 &  bit2 &  bit3 & ~bit4), # 7: 1st, 2nd, 3rd bits
              :(~bit1 & ~bit2 & ~bit3 &  bit4)) # 8: 4th bit

    # Generator with the on bits of cells that should stay alive
    alives = (onbit for (i, onbit) in enumerate(onbits) if i in rulesums(N))

    # Return the appropriate pieces OR'ed together, or zero(T) if the rule is "always die"
    return isempty(alives) ? :(zero(T)) : Expr(:call, :|, alives...)
end



"""
    updatedcluster(above::Integer, current::Integer, below::Integer, ::LifeRule{B, S})

Update and return `current` given adjacent clusters `above` and `below` using bit-twiddling.

`updatedcluster` uses bitwise operations to efficiently count the number of neighbors of
each cell in the `current` cluster and return the updated cluster. It employs half and full
adders to calculate the sums and carries needed.

Updates are done according to the specified [`LifeRule`](@ref), which contains birth (`B`)
and survival (`S`) conditions.
"""
function updatedcluster(above, current, below, ::LifeRule{B, S}) where {B, S}
    bit1, bit2, bit3, bit4 = bitsums(above, current, below)

    # Update current according to the survival and birth rules
    current &= alivebits(bit1, bit2, bit3, bit4, S) # survival
    current |= alivebits(bit1, bit2, bit3, bit4, B) # birth

    # Return the updated cluster
    return current
end



#=
# An example of specializing updatedcluster for a particular rule (B3/S23)
# This doesn't seem to be needed, the compiler does great at optimizing the generic version
function updatedcluster(above, current, below, ::LifeRule{Rule(3), Rule(2, 3)})
    bit1, bit2, bit3, bit4 = bitsums(above, current, below)

    current |=  bit1 # birth 1, 3, 5, and 7
    current &=  bit2 # kill 1, 4, 5, and 8
    current &= ~bit3 # kill 6 and 7

    return current
end

# updatedcluster specialization for B36/S23 (high life)
function updatedcluster(above, current, below, ::LifeRule{Rule(3, 6), Rule(2, 3)})
    bit1, bit2, bit3, bit4 = bitsums(above, current, below)

    current &=  bit2 # kill 1, 4, 5, and 8
    current &= ~bit3 # kill 6 and 7
    current |=  bit2 & (bit1 ⊻ bit3) # birth 3 or 6, but not 7

    return current
end
=#
