#!/bin/bash

#=
exec julia --threads=4 --project="$(dirname "$(realpath "$0")")/.." "$0"
=#

using LifeGame, Plots, BenchmarkTools



# Return step! cell updates per nanosecond for each side length in sidelengths
function benchmarkstep(sidelengths, par)
    # Warm up the CPU for a minute
    t = time()
    while time()-t < 60
        step!(LifeGrid(10_000, 10_000), par ? parallel : serial)
    end

    return [
        begin
            # Free up memory
            GC.gc()

            # Create the grid, allowing each core some work if par
            lg = if par
                w = N/Threads.nthreads()
                CType = w <= 6 ? UInt8 : w <= 14 ? UInt16 : w <= 30 ? UInt32 : UInt64
                LifeGrid(N, N, CType = CType)
            else
                LifeGrid(N, N)
            end

            # Get benchmark results
            result = @benchmark step!($lg, $par ? parallel : parallel)

            # Return cells per ns
            N^2/mean(result.times)
        end for N in sidelengths
    ]
end



# Get dense benchmarks of small sizes, sparse benchmarks of bigger sizes
sidelengths = (small = 1:500, large = 5 .* 2 .^ (1:16))
results = (
    smallserial = benchmarkstep(sidelengths.small, false),
    smallparallel = benchmarkstep(sidelengths.small, true),
    largeserial = benchmarkstep(sidelengths.large, false),
    largeparallel = benchmarkstep(sidelengths.large, true),
)



# Plot and save results
smallplotkw = (
    title = "Cells updated per nanosecond",
    legend_position = :topleft,
    xlabel = "LifeGrid side length",
    ylabel = "Updates/ns",
    margin = (5, :mm),
    size = (600, 400),
)
largeplotkw = (
    xscale = :log10,
    xticks = (sidelengths.large, sidelengths.large),
    xrotation = 45,
    marker = :circle,
    markerstrokewidth = 0,
    smallplotkw...,
)

smallplot = plot(sidelengths.small, results.smallserial; label = "serial", smallplotkw...)
plot!(
    smallplot,
    sidelengths.small,
    results.smallparallel;
    label = "parallel",
    smallplotkw...,
)
savefig(smallplot, "benchmark-results-small.png")

largeplot = plot(sidelengths.large, results.largeserial; label = "serial", largeplotkw...)
plot!(
    largeplot,
    sidelengths.large,
    results.largeparallel;
    label = "parallel",
    largeplotkw...,
)
savefig(largeplot, "benchmark-results-large.png")
