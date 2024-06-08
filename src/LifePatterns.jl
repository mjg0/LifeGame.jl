export LifePatterns

"""
    LifePatterns

A collection of common [`LifePattern`](@ref)s that can be inserted into [`LifeGrid`](@ref)s.
"""
module LifePatterns



import ..LifeGame.LifePattern



# Macro to tersely create a constant pattern
macro pattern(name, pattern)
    return quote
        const $(esc(name)) = LifePattern($pattern)
    end
end



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