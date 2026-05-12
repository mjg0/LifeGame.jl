using Test, Crayons, Random, LifeGame



include("utils.jl")

include("SlowLifeGrid.jl")

@testset "LifeGame" begin
    #include("constructor.jl")
    #include("indexing.jl")
    #include("updatedcluster.jl")
    #include("step.jl")
    #include("io.jl")
    include("patterns.jl")
end