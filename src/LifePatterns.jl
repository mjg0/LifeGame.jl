export LifePatterns



"""
    LifePatterns

A collection of common patterns meant for insertion into [`LifeGrid`](@ref)s.
"""
module LifePatterns



using ..LifeGame



# Macro to concisely create a constant pattern
macro pattern(name, pattern)
    return quote
        const $(esc(name)) = LifePattern($pattern.!=0)
    end
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