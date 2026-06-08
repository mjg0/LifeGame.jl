using Test, Crayons, Random, LifeGame

include("SlowLifeGrid.jl")



@testset "LifeGame" begin
    include("utils.jl")

    include("updatedcluster.jl")

    include("constructors.jl")

    include("mutators.jl")

    include("io.jl")
end
