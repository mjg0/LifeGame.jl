begin
    rng = MersenneTwister(1)

    N = 100
    iters = 100

    # Ensure we get the rules created by @specialize_updatedcluster
    specialized = (LifeGrid(rand(rng, Bool, rand(rng, 10:200), rand(rng, 10:200)); rule=R)
                   for R in ("B3/S23", "B36/S23", "B2/S", "B234/S"))

    grids = (specialized..., ntuple(_ -> randlifegrid(rng), N - length(specialized))...)
    references = ntuple(i -> SlowLifeGrid(grids[i]; rule=repr(rule(grids[i]))), N)

    @testset "step!" begin
        for (grid, reference) in zip(grids, references), _ in 1:iters
            stepandcheck!(reference, grid, rng)
        end
    end

    @testset "setindex!" begin
        for (grid, reference) in zip(grids, references), _ in 1:iters
            indices = if rand(rng, Bool) # dense
                rand(rng, CartesianIndices(grid), cld(length(grid), 4))
            else # sparse
                [rand(rng, CartesianIndices(grid))]
            end
            mutate!(A) = for I in indices
                A[I] = !A[I]
            end
            I = rand(rng, CartesianIndices(grid))
            stepandcheck!(reference, grid, rng, mutate!)
        end
    end

    @testset "insert!" begin
        for (grid, reference) in zip(grids, references), _ in 1:iters
            M, N = size(grid)
            p = rand(rng, Bool, rand(rng, 1:M), rand(rng, 1:N))
            
            I = rand(rng, CartesianIndex((1, 1)):CartesianIndex(size(reference).-size(p).+1))
            mutate!(A) = insert!(A, I, rand(rng, Bool) ? p : LifePattern(p))
            stepandcheck!(reference, grid, rng, mutate!)
        end
    end
end