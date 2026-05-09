export LifePatterns



"""
    LifePatterns

A collection of common patterns meant for insertion into [`LifeGrid`](@ref)s.
"""
module LifePatterns



using ..LifeGame



"""
    @pattern name pattern

Concisely create constant `LifePattern`s.

The name of the resultant constant(s) will be `name` if the supplied `pattern` array is
invariant under rotation and reflection, or `name<N>` otherwise. Each unique
rotation/reflection combination on `pattern` will correspond to a different `<N>`, where
`N∈[1,8]`.

# Examples

This will result in two constants, `LifePattern.H1` and `LifePattern.H2`:

```julia
@pattern H [1 0 1
            1 1 1
            1 0 1]
```
"""
macro pattern(name, pattern)
    # Turn pattern into a literal matrix
    P = Bool.(Core.eval(__module__, :($pattern .!= 0)))

    # Turn a literal matrix back into an Expr
    matrixexpr(x) = Expr(:vcat, [
        Expr(:row, [x[i, j] for j in axes(x, 2)]...)
        for i in axes(x, 1)
    ]...)

    # Find all unique rotation/reflection combos
    patterns = Matrix{Bool}[]
    for rotation in (x->x, rotl90, rot180, rotr90)
        for reflection in (x->x, x->reverse(x; dims=2))
            transformed = reflection(rotation(P))
            if transformed ∉ patterns
                push!(patterns, transformed)
            end
        end
    end

    # Create and return the list of consts
    defs = if length(patterns) == 1
        [:(const $(Symbol(name))    = LifePattern($(matrixexpr(first(patterns)));
                                                  mutable=false))]
    else
        [:(const $(Symbol(name, k)) = LifePattern($(matrixexpr(p));
                                                  mutable=false))
         for (k, p) in enumerate(patterns)]
    end
    return esc(Expr(:block, defs...))
end



# Patterns
@pattern block          [1 1
                         1 1]

@pattern beehive        [0 1 1 0
                         1 0 0 1
                         0 1 1 0]

@pattern blinker        [0 1 0
                         0 1 0
                         0 1 0]

@pattern toad           [0 0 0 0
                         0 1 1 1
                         1 1 1 0
                         0 0 0 0]

@pattern beacon         [1 1 0 0
                         1 1 0 0
                         0 0 1 1
                         0 0 1 1]

@pattern pulsar         [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                         0 0 0 1 1 1 0 0 0 1 1 1 0 0 0
                         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                         0 1 0 0 0 0 1 0 1 0 0 0 0 1 0
                         0 1 0 0 0 0 1 0 1 0 0 0 0 1 0
                         0 1 0 0 0 0 1 0 1 0 0 0 0 1 0
                         0 0 0 1 1 1 0 0 0 1 1 1 0 0 0
                         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                         0 0 0 1 1 1 0 0 0 1 1 1 0 0 0
                         0 1 0 0 0 0 1 0 1 0 0 0 0 1 0
                         0 1 0 0 0 0 1 0 1 0 0 0 0 1 0
                         0 1 0 0 0 0 1 0 1 0 0 0 0 1 0
                         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                         0 0 0 1 1 1 0 0 0 1 1 1 0 0 0
                         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]

@pattern pentadecathlon [0 0 0 0 0 0 0 0 0
                         0 0 0 0 0 0 0 0 0
                         0 0 0 0 0 0 0 0 0
                         0 0 0 0 1 0 0 0 0
                         0 0 0 0 1 0 0 0 0
                         0 0 0 1 0 1 0 0 0
                         0 0 0 0 1 0 0 0 0
                         0 0 0 0 1 0 0 0 0
                         0 0 0 0 1 0 0 0 0
                         0 0 0 0 1 0 0 0 0
                         0 0 0 1 0 1 0 0 0
                         0 0 0 0 1 0 0 0 0
                         0 0 0 0 1 0 0 0 0
                         0 0 0 0 0 0 0 0 0
                         0 0 0 0 0 0 0 0 0
                         0 0 0 0 0 0 0 0 0]

@pattern glider         [0 1 0
                         0 0 1
                         1 1 1]



end # module