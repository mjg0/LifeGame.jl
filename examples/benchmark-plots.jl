#!/bin/bash

#=
exec julia --threads=4 --project="$(dirname "$(realpath "$0")")/.." "$0"
=#

using LifeGame, Plots, BenchmarkTools, Random, Statistics



# @benchmark is VERY picky about being called just right
@generated function meansteptimesns(lg, args::Vararg{Type, N}) where N
    modes = map(T -> T.parameters[1], args)

    return quote
        mean(@benchmark step!($(Expr(:$, :lg)), $(modes...))).time
    end
end



# Return step! cell updates per nanosecond for each side length in sidelengths
function benchmarkstep(sidelengths, threadmode, sparsity)
    par = threadmode == :serial ? serial : parallel
    # Warm up the CPU for a minute
    t = time()
    while time() - t < 10
        step!(LifeGrid(10_000, 10_000), par)
    end

    return [
        begin
            # Free up memory
            GC.gc()

            # Create the LifeGrid
            lg = if threadmode == :parallel
                w = N / Threads.nthreads()
                C = w <= 6 ? UInt8 : w <= 14 ? UInt16 : w <= 30 ? UInt32 : UInt64
                LifeGrid(N, N; CType=C)
            else
                LifeGrid(N, N)
            end
            rand!(lg) # Otherwise the sparse algorithm doesn't have to do any work

            # Get benchmark results
            meantime = if threadmode == :default # don't specify parallel or serial
                meansteptimesns(lg, sparsity)
            else
                meansteptimesns(lg, sparsity, par)
            end

            # Return cells per ns
            N^2 / meantime
        end for N in sidelengths
    ]
end



# Get dense benchmarks of small sizes, plus both sparse and dense for bigger sizes
sidelengths = (small = 1:400, large = 15 .* 2 .^ (1:14))
defaultlengths = (small = 1:150, large = 15 .* 2 .^ (1:5))

series = (
    (key = :serial, label = "serial", threadmode = :serial, lengths = sidelengths),
    (key = :parallel, label = "parallel", threadmode = :parallel, lengths = sidelengths),
    (key = :default, label = "default", threadmode = :default, lengths = defaultlengths),
)

benchmarkset(which, sparsity) = (
    serial = benchmarkstep(series[1].lengths[which], series[1].threadmode, sparsity),
    parallel = benchmarkstep(series[2].lengths[which], series[2].threadmode, sparsity),
    default = benchmarkstep(series[3].lengths[which], series[3].threadmode, sparsity),
)

results = (
    small = benchmarkset(:small, dense),
    large_sparse = benchmarkset(:large, sparse),
    large_dense = benchmarkset(:large, dense),
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

function addplotseries!(p, resultset, which, plotkw)
    for entry in series
        plot!(
            p,
            entry.lengths[which],
            resultset[entry.key];
            label = entry.label,
            plotkw...,
        )
    end
    p
end

function overlay_sparse!(p)
    for (i, entry) in pairs(series)
        plot!(
            p,
            entry.lengths.large,
            results.large_sparse[entry.key];
            label = "",
            color = i,
            linewidth = 0.5,
            linestyle = :dash,
            marker = :none,
        )
    end
    p
end

smallplot = plot(; smallplotkw...)
addplotseries!(smallplot, results.small, :small, smallplotkw)
savefig(smallplot, "benchmark-results-small.png")

largeplot = plot(; largeplotkw...)
addplotseries!(largeplot, results.large_dense, :large, largeplotkw)
overlay_sparse!(largeplot)
savefig(largeplot, "benchmark-results-large.png")
