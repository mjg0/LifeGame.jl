export LifePatterns



# Optimization of insert! for sparse matrices
Base.@propagate_inbounds function
Base.insert!(lg::LifeGrid, I::CartesianIndex{2}, pattern::SparseMatrixCSC)
    rows = rowvals(pattern)

    for j in axes(pattern, 2)
        for nz in nzrange(pattern, j)
            i = rows[nz]
            lg[I+CartesianIndex(i, j)-one(I)] = true
        end
    end

    return lg
end

# Other matrices use `findall`, which is slow and allocating
Base.@propagate_inbounds function
Base.insert!(lg::LifeGrid, I::CartesianIndex{2}, pattern::AbstractMatrix{<:Number})
    for i in CartesianIndices(pattern)
        if pattern[i] != zero(eltype(pattern))
            lg[I+i-one(I)] = true
        end
    end

    return lg
end

Base.@propagate_inbounds function Base.insert!(lg::LifeGrid, i::Integer, j::Integer, pattern)
    return insert!(lg, CartesianIndex(i, j), pattern)
end



"""
    LifePatterns

A collection of common patterns meant for insertion into [`LifeGrid`](@ref)s.
"""
module LifePatterns



using SparseArrays



# Macro to concisely create a constant pattern
macro pattern(name, pattern)
    return quote
        const $(esc(name)) = sparse($pattern.!=0)
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