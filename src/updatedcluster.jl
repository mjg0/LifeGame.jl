"""
    bitsums(above, current, below)

Return a 4-tuple containing the bitwise sum of the neighbor bits of the cluster `current`.

The first element of the returned tuple represents the ones place, the second the twos
place, the third the fours place, and the fourth the eights place.

# Extended help

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

# Extended help

It's possible to specialize `updatedcluster` for a specific rule. The compiler does a good
job with the generic `updatedcluster` so usually the performance gain is marginal, but for
some rules a modest boost (~10%) is attainable. The signature for the specialization should
be:

```julia
LifeGame.updatedcluster(above, current, below,
                        ::LifeGame.LifeRule{LifeGame.Rule(B...), LifeGame.Rule(S...)})
```

...where `B...` and `S...` are lists of neighbor sums that should result in birth and
survival, respectively. For example, to specialize `updatedcluster` for `B123/S45`, the
signature is:

```julia
LifeGame.updatedcluster(above, current, below,
                        ::LifeGame.LifeRule{LifeGame.Rule(1, 2, 3), LifeGame.Rule(4, 5)})
```

`above`, `current`, and `below` are "clusters": 64-bit unsigned integers representing rows
of 62 living cells, with a single cell of padding at each end to allow single bitshifts not
to roll relevant cells off of the edges. `above` is the cluster above `current`, and `below`
is the cluster below. `updatedcluster` returns a single 64-bit unsingned integer containing
the 62 cells (the first and last bits are ignored) of `current` stepped forward one
generation.

`LifeGame.bitsums` will be helpful for most specializations. It takes `above`, `current`,
and `below` and returns the bitwise sum as 4 64-bit unsigned integers: the first contains
the ones place, the second the twos place, the third the fours place, and the fourth the
eights place.

As an example, the specialization for the rule `B3/S4` might look like:

```julia
function LifeGame.updatedcluster(above, current, below,
                                 ::LifeGame.LifeRule{LifeGame.Rule(1), LifeGame.Rule(2)})
    # Get the bitwise sum of the neighbor cells of each bit in `current`
    1s_place, 2s_pace, 4s_place, 8s_place = LifeGame.bitsums(above, current, below)

    # Kill all cells but those with a neighbor sum of 4
    current &= ~1s_place & ~2s_place & 4s_place & ~8s_place

    # Birth cells that had 3 neighbors
    current |= 1s_place & 2s_place & ~4s_place & ~8s_place

    # Return the updated cluster
    return current
end
```
"""
function updatedcluster(above, current, below, ::LifeRule{B, S}) where {B, S}
    bit1, bit2, bit3, bit4 = bitsums(above, current, below)

    # Update current according to the survival and birth rules
    current &= alivebits(bit1, bit2, bit3, bit4, S) # survival
    current |= alivebits(bit1, bit2, bit3, bit4, B) # birth

    # Return the updated cluster
    return current
end



"""
    @specialize_updatedcluster birth survival body

Create a specialization of [`updatedcluster`](@ref) for the given rule.

`birth` and `survival` are lists of neighbor sums causing birth and survival, respectively.

`body` is a statement that can use `above`, `below`, `current` (the arguments to
`updatedcluster`), `bit1`, `bit2`, `bit3`, and `bit4` (the results of
`bitsums(above, below, current)`) to compute the updated state of `current`.

As an example, the `updatedcluster` specialization for `B2/S` is created with:

```julia
@specialize_updatedcluster [2] [] ~bit1 & bit2 & ~bit3
```
"""
macro specialize_updatedcluster(birth, survival, body)
    liferule = :(::LifeRule{Rule($(birth)...), Rule($(survival)...)})
    qualified_name = esc(GlobalRef(LifeGame, :updatedcluster))

    return quote
        function $qualified_name(above, current, below, $liferule)
            bit1, bit2, bit3, bit4 = bitsums(above, current, below)
            $body
        end
    end
end



# Specializations
# B3/S23 (Conway's life)
@specialize_updatedcluster [3]       [2, 3] (current | bit1) & bit2 & ~bit3

# B36/S23 (highlife)
@specialize_updatedcluster [3, 6]    [2, 3] current & bit2 & ~bit3 | bit2 & (bit1 ⊻ bit3)

# B2/S (seeds)
@specialize_updatedcluster [2]       []     ~bit1 & bit2 & ~bit3

# B234/S (Persian rug)
@specialize_updatedcluster [2, 3, 4] []     bit2 & ~bit3 | ~bit1 & ~bit2 & bit3
