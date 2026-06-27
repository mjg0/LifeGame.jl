export step!, dense, sparse, serial, parallel



# @batch is noticeably expensive for small grids, this allows that to be avoided
abstract type Threaded end
struct serial <: Threaded end
struct parallel <: Threaded end



# Type to indicate whether the dense or sparse algorithm should be used
abstract type Sparsity end
struct dense <: Sparsity end
struct sparse <: Sparsity end



# Check whether a rule causes birth with zero surviving neighbors
function haszerobirth(::Type{LifeRule{B,S}}) where {B,S}
    return B & 0x0001 == 0x0001
end



"""
    step!(lg::LifeGrid, sparseordense, serialorparallel)
    step!(lg::LifeGrid, serialorparallel, sparseordense)

Update `lg` one generation according to the [`rule`](@ref) associated with it.

All cells outside of the grid boundary are fixed at zero.

`sparseordense` can be `sparse` or `dense`; `serialorparallel` can be `serial` or
`parallel`. They default to `dense` and `parallel` respectively.
"""
@inline function step!(
    lg::LifeGrid{R,C,H,Tall,Wide},
    ::Type{Sp},
    ::Type{Th},
) where {R,C,H,Tall,Wide,Sp<:Sparsity,Th<:Threaded}
    # If B0 is part of the rule, the sparse algorithm can't be used
    sparseordense = Sp === sparse && haszerobirth(R) ? dense : Sp

    # Only parallelize if there's a need
    if Wide && Th === parallel
        # Manual thread bounds so buffers don't get mixed
        nthreads = min(Threads.nthreads(), size(lg.grid, 2))
        colsperthread = cld(size(lg.grid, 2), nthreads)

        # Outer loop over threads, inner over columns
        @inbounds @batch for tid = 1:nthreads
            J1 = colsperthread * (tid - 1) + 1
            J2 = min(colsperthread * tid, size(lg.grid, 2))
            for J = J1:J2
                updatecolumn!(lg, J, tid, sparseordense)
            end
        end
    else
        @inbounds for J in axes(lg.grid, 2)
            updatecolumn!(lg, J, 1, sparseordense)
        end
    end

    # Swap halos
    lg.halos.currentleft, lg.halos.nextleft = lg.halos.nextleft, lg.halos.currentleft
    lg.halos.currentright, lg.halos.nextright = lg.halos.nextright, lg.halos.currentright

    # Update sparse activity trackers
    if sparseordense == sparse
        lg.changed.current, lg.changed.next = lg.changed.next, lg.changed.current
        lg.changed.allcurrent = false
    else
        lg.changed.allcurrent = true
    end

    return lg
end

step!(lg::LifeGrid, ::Type{Sp}) where {Sp<:Sparsity} = step!(lg, Sp, parallel)
step!(lg::LifeGrid, ::Type{Th}) where {Th<:Threaded} = step!(lg, dense, Th)
step!(lg::LifeGrid, ::Type{Th}, ::Type{Sp}) where {Th<:Threaded,Sp<:Sparsity} =
    step!(lg, Sp, Th)
step!(lg::LifeGrid) = step!(lg, parallel)
