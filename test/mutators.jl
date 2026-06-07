begin
    rng = MersenneTwister(1)

    N = 100

    grids = ntuple(_ -> randlifegrid(rng), N)
    references = ntuple(i -> SlowLifeGrid(grids[i]; rule=repr(rule(grids[i]))), N)

    @testset "step!" begin
        for (grid, reference) in zip(grids, references), _ in 1:5
            stepandcheck!(reference, grid, rng)
        end
    end

    @testset "setindex!" begin
        for (grid, reference) in zip(grids, references), _ in 1:5
            indices = rand(rng, CartesianIndices(grid), cld(length(grid), 4))
            mutate!(A) = for I in indices
                A[I] = !A[I]
            end
            I = rand(rng, CartesianIndices(grid))
            stepandcheck!(reference, grid, rng, mutate!)
        end
    end

    @testset "insert!" begin
        for (grid, reference) in zip(grids, references), _ in 1:5
            M, N = size(grid)
            p = rand(rng, Bool, rand(rng, 1:M), rand(rng, 1:N))
            
            I = rand(rng, CartesianIndex((1, 1)):CartesianIndex(size(reference).-size(p).+1))
            mutate!(A) = insert!(A, I, rand(rng, Bool) ? p : LifePattern(p))
            stepandcheck!(reference, grid, rng, mutate!)
        end
    end
end