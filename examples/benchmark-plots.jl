#!/bin/bash

#=
exec julia --threads=4 --project="$(dirname "$(realpath "$0")")/.." "$0"
=#

using LifeGame, Plots, BenchmarkTools



# Return step! cell updates per nanosecond for each side length in sidelengths
function benchmarkstep(sidelengths, type)
    # Warm up the CPU for a minute
    t = time()
    while time() - t < 10
        step!(LifeGrid(10_000, 10_000), type == :serial ? serial : parallel)
    end

    return [
        begin
            # Free up memory
            GC.gc()

            # Get benchmark results
            result = if type == :parallel
                w = N / Threads.nthreads()
                CType = w <= 6 ? UInt8 : w <= 14 ? UInt16 : w <= 30 ? UInt32 : UInt64
                lg = LifeGrid(N, N, CType=CType)
                @benchmark step!($lg)
            elseif type == :serial
                lg = LifeGrid(N, N)
                @benchmark step!($lg, serial)
            else # :default
                lg = LifeGrid(N, N)
                @benchmark step!($lg)
            end

            # Return cells per ns
            N^2 / mean(result.times)
        end for N in sidelengths
    ]
end



# Get dense benchmarks of small sizes, sparse benchmarks of bigger sizes
sidelengths = (small=1:500, large=5 .* 2 .^ (1:16))
defaultlengths = (small=1:150, large=5 .* 2 .^ (1:6))
results = (
    smallserial=benchmarkstep(sidelengths.small, :serial),
    smallparallel=benchmarkstep(sidelengths.small, :parallel),
    smalldefault=benchmarkstep(defaultlengths.small, :default),
    largeserial=benchmarkstep(sidelengths.large, :serial),
    largeparallel=benchmarkstep(sidelengths.large, :parallel),
    largedefault=benchmarkstep(defaultlengths.large, :default),
)



# Plot and save results
smallplotkw = (
    title="Cells updated per nanosecond",
    legend_position=:topleft,
    xlabel="LifeGrid side length",
    ylabel="Updates/ns",
    margin=(5, :mm),
    size=(600, 400),
)
largeplotkw = (
    xscale=:log10,
    xticks=(sidelengths.large, sidelengths.large),
    xrotation=45,
    marker=:circle,
    markerstrokewidth=0,
    smallplotkw...,
)

smallplot = plot(sidelengths.small, results.smallserial; label="serial", smallplotkw...)
plot!(smallplot, sidelengths.small, results.smallparallel; label="parallel", smallplotkw...)
plot!(smallplot, defaultlengths.small, results.smalldefault; label="default", smallplotkw...)
savefig(smallplot, "benchmark-results-small.png")

largeplot = plot(sidelengths.large, results.largeserial; label="serial", largeplotkw...)
plot!(largeplot, sidelengths.large, results.largeparallel; label="parallel", largeplotkw...)
plot!(largeplot, defaultlengths.large, results.largedefault; label="default", largeplotkw...)
savefig(largeplot, "benchmark-results-large.png")
